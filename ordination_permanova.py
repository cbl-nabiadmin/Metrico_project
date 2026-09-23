#!/usr/bin/env python3
"""
Ordination and PERMANOVA for the Metrico Sample Set
--------------------------------------------------
Reads a wide feature matrix (KO, EC, COG, CAZy) and metadata to calculate:
1. PCoA Ordination (Aitchison CLR or Bray-Curtis)
2. PERMANOVA for host, batch, and stratum-restricted host effects
3. Static PNG and Interactive HTML scatter plots
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd

# Configure Matplotlib for non-interactive rendering in headless Docker environments
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.express as px

# Safe fallback import for metrico_checks module
try:
    import metrico_checks as mc
    HAS_MC = True
except ImportError:
    HAS_MC = False

# =========================================================
# 1. COMMAND LINE ARGUMENTS & PATH RESOLUTION
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")

parser = argparse.ArgumentParser(
    description="PCoA ordination and PERMANOVA for Metrico sample matrices",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-i", "--input", type=str, default="all_ko_matrix.tsv", help="Wide feature matrix file")
parser.add_argument("-m", "--metadata", type=str, default="Plant-host.csv", help="Metadata CSV file")
parser.add_argument("-o", "--output-dir", type=str, default="/Metrico_project/output", help="Output directory")
parser.add_argument("--prefix", type=str, default="Ordination", help="Prefix for output files")
parser.add_argument("--group", type=str, default="Host", help="Metadata grouping column")
parser.add_argument("--batch", type=str, default="Study", help="Metadata batch column")
parser.add_argument("--metric", choices=["aitchison", "bray"], default="aitchison", help="Distance metric")
parser.add_argument("--permutations", type=int, default=999, help="PERMANOVA permutations")
parser.add_argument("--min-group", type=int, default=3, help="Minimum sample size per group")
parser.add_argument("--min-prevalence", type=float, default=0.1, help="Minimum feature prevalence fraction")
parser.add_argument("--ellipse-level", type=float, default=0.95, help="Group confidence ellipse coverage")
parser.add_argument("--ellipse-min-n", type=int, default=6, help="Minimum sample size for drawing ellipses")
parser.add_argument("--seed", type=int, default=1234, help="Random seed")

args = parser.parse_args()

OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)
os.makedirs(OUTPUT_DIR, exist_ok=True)
rng = np.random.default_rng(args.seed)

matrix_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
meta_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

for path, description in [(matrix_file, "Feature matrix"), (meta_file, "Metadata")]:
    if not os.path.exists(path):
        print(f"❌ ERROR: {description} not found at {path}", file=sys.stderr)
        sys.exit(1)

_base = os.path.basename(matrix_file).lower()
if "ko" in _base:        FEATURE_LABEL = "KO IDs"
elif "ec" in _base:      FEATURE_LABEL = "EC Numbers"
elif "cog" in _base:     FEATURE_LABEL = "COG Categories"
elif "cazy" in _base:    FEATURE_LABEL = "CAZy Classes"
else:                    FEATURE_LABEL = "Functional Profiles"

print(f"📂 Matrix File   : {matrix_file}")
print(f"📂 Metadata File : {meta_file}")
print(f"📊 Distance Metric: {args.metric}")

# =========================================================
# 2. LOAD MATRIX AND METADATA
# =========================================================
try:
    df = pd.read_csv(matrix_file, sep=None, engine="python", index_col=0)
except Exception:
    df = pd.read_csv(matrix_file, sep="\t", index_col=0)

df = df.loc[~df.index.astype(str).isin(["EC", "KO", "Sample", "SampleID"])]
df = df.loc[:, ~df.columns.astype(str).isin(["EC", "KO", "Sample", "SampleID"])]
df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

meta = pd.read_csv(meta_file)
meta.columns = meta.columns.str.strip()
s_col = next((c for c in meta.columns if "sample" in c.lower()), meta.columns[0])
g_col = args.group if args.group in meta.columns else next((c for c in meta.columns if "host" in c.lower()), None)

if g_col is None:
    print(f"❌ ERROR: Grouping column '{args.group}' not found in metadata", file=sys.stderr)
    sys.exit(1)

meta[s_col] = meta[s_col].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta = meta[meta[s_col].str.lower() != "sampleid"].copy()
meta = meta.drop_duplicates(subset=[s_col]).set_index(s_col)

ids = set(meta.index)
cols_as_samples = len(ids.intersection(set(df.columns.astype(str))))
rows_as_samples = len(ids.intersection(set(df.index.astype(str))))

if cols_as_samples >= rows_as_samples:
    df = df.T
    print(f"🔄 Matrix transposed: samples were columns ({cols_as_samples} matched)")
else:
    print(f"📊 Matrix aligned as samples x features ({rows_as_samples} matched)")

df.index = df.index.astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)

common = sorted(set(df.index).intersection(ids))
if HAS_MC:
    mc.check_overlap(len(common), meta.index.nunique(), "samples", data_total=df.index.nunique())

df = df.loc[common]
meta = meta.loc[common]

# =========================================================
# 3. FILTERING
# =========================================================
counts = meta[g_col].value_counts()
keep_groups = counts[counts >= args.min_group].index
dropped = counts[counts < args.min_group]

if len(dropped):
    print(f"⚠️ Dropping {len(dropped)} group(s) with n < {args.min_group}: " + ", ".join(f"{g} (n={n})" for g, n in dropped.items()))

mask = meta[g_col].isin(keep_groups)
df, meta = df.loc[mask.values], meta.loc[mask.values]

prev = (df > 0).mean(axis=0)
df = df.loc[:, prev >= args.min_prevalence]
print(f"📊 Analyzing {df.shape[0]} samples x {df.shape[1]} features across {len(keep_groups)} '{g_col}' groups")

if df.shape[0] < 6:
    print("❌ ERROR: Too few samples remaining for ordination analysis.", file=sys.stderr)
    sys.exit(1)

# =========================================================
# 4. DISTANCE TRANSFORMATION
# =========================================================
X = df.to_numpy(dtype=float)

if args.metric == "aitchison":
    X = X + 0.5
    X = X / X.sum(axis=1, keepdims=True)
    logX = np.log(X)
    clr = logX - logX.mean(axis=1, keepdims=True)
    diff = clr[:, None, :] - clr[None, :, :]
    D = np.sqrt((diff ** 2).sum(axis=-1))
    metric_label = "Aitchison distance (Euclidean on CLR)"
else:
    rel = X / np.where(X.sum(axis=1, keepdims=True) == 0, 1, X.sum(axis=1, keepdims=True))
    s = rel[:, None, :] + rel[None, :, :]
    d = np.abs(rel[:, None, :] - rel[None, :, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        D = np.where(s.sum(axis=-1) == 0, 0.0, d.sum(axis=-1) / s.sum(axis=-1))
    metric_label = "Bray-Curtis distance"

np.fill_diagonal(D, 0.0)

# =========================================================
# 5. PCoA CALCULATION
# =========================================================
def pcoa(D, n_axes=3):
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    B = (B + B.T) / 2.0
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals > 0
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    explained = vals[pos] / vals[pos].sum()
    return coords[:, :n_axes], explained[:n_axes], explained

coords, expl, expl_all = pcoa(D)
print(f"📊 PCoA Axes Variance Explained: PCo1 ({expl[0]:.1%}), PCo2 ({expl[1]:.1%}), PCo3 ({expl[2]:.1%})")

ord_df = pd.DataFrame(coords, columns=["PCo1", "PCo2", "PCo3"], index=df.index)
ord_df[g_col] = meta[g_col].values
if args.batch in meta.columns:
    ord_df[args.batch] = meta[args.batch].values

# =========================================================
# 6. PERMANOVA CALCULATIONS
# =========================================================
def permanova(D, labels, permutations=999, strata=None, rng=None):
    rng = rng or np.random.default_rng(0)
    labels = np.asarray(labels)
    n = D.shape[0]
    D2 = D ** 2
    groups = np.unique(labels)
    a = len(groups)
    if a < 2:
        return dict(error="Fewer than two groups")

    def ss(labels_):
        ss_total = D2.sum() / (2.0 * n)
        ss_within = 0.0
        for g in np.unique(labels_):
            idx = np.where(labels_ == g)[0]
            ng = len(idx)
            if ng > 1:
                ss_within += D2[np.ix_(idx, idx)].sum() / (2.0 * ng)
        return ss_total, ss_within

    ss_total, ss_within = ss(labels)
    ss_between = ss_total - ss_within
    denom = ss_within / (n - a) if n > a else np.nan
    F = (ss_between / (a - 1)) / denom if denom and denom > 0 else np.nan
    R2 = ss_between / ss_total if ss_total > 0 else np.nan

    if strata is None:
        perm_idx = [rng.permutation(n) for _ in range(permutations)]
    else:
        strata = np.asarray(strata)
        perm_idx = []
        for _ in range(permutations):
            p = np.arange(n)
            for s in np.unique(strata):
                pos = np.where(strata == s)[0]
                p[pos] = rng.permutation(pos)
            perm_idx.append(p)

    count = 0
    for p in perm_idx:
        _, sw = ss(labels[p])
        sb = ss_total - sw
        dn = sw / (n - a)
        Fp = (sb / (a - 1)) / dn if dn > 0 else np.nan
        if np.isfinite(Fp) and Fp >= F:
            count += 1
    pval = (count + 1) / (permutations + 1)

    return dict(groups=int(a), n=int(n), pseudo_F=float(F), R2=float(R2),
                p_value=float(pval), permutations=int(permutations),
                restricted_within_strata=strata is not None)

results = {}
results[g_col] = permanova(D, meta[g_col].values, args.permutations, rng=rng)
r = results[g_col]
print(f"📊 PERMANOVA ({g_col}): R2={r['R2']:.3f}, pseudo-F={r['pseudo_F']:.2f}, p={r['p_value']:.4f}")

if args.batch in meta.columns:
    _unres = meta[args.batch].isna() | meta[args.batch].astype(str).str.strip().isin(["", "UNKNOWN", "unknown", "NA", "nan", "None"])
    _n_unres = int(_unres.sum())
    if _n_unres:
        meta.loc[_unres, args.batch] = ["UNRESOLVED_%s" % sid for sid in meta.index[_unres]]
        if HAS_MC:
            mc.record(unresolved_batch_samples=_n_unres)

has_batch = args.batch in meta.columns and meta[args.batch].nunique() > 1
if has_batch:
    results[args.batch] = permanova(D, meta[args.batch].values, args.permutations, rng=rng)
    results[f"{g_col}_within_{args.batch}"] = permanova(D, meta[g_col].values, args.permutations, strata=meta[args.batch].values, rng=rng)

# =========================================================
# 7. PLOTS GENERATION
# =========================================================
def ellipse_points(x, y, level=0.95, n_points=120, min_n=3):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < max(3, min_n):
        return None
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)):
        return None
    vals, vecs = np.linalg.eigh(cov)
    if vals.min() <= 1e-12:
        return None
    radius = np.sqrt(-2.0 * np.log(1.0 - level))
    t = np.linspace(0.0, 2.0 * np.pi, n_points)
    circle = np.stack([np.cos(t), np.sin(t)])
    ell = vecs @ np.diag(np.sqrt(vals) * radius) @ circle
    return ell[0] + x.mean(), ell[1] + y.mean()

def scatter_png(color_col, fname, title, subtitle_text=""):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    levels = sorted(ord_df[color_col].astype(str).unique())
    cmap = plt.get_cmap("tab20" if len(levels) > 10 else "tab10")
    
    for i, lev in enumerate(levels):
        sel = ord_df[ord_df[color_col].astype(str) == lev]
        colour = cmap(i % cmap.N)

        if args.ellipse_level > 0:
            ell = ellipse_points(sel["PCo1"], sel["PCo2"], args.ellipse_level, min_n=args.ellipse_min_n)
            if ell is not None:
                ax.fill(ell[0], ell[1], color=colour, alpha=0.11, linewidth=0, zorder=1)
                ax.plot(ell[0], ell[1], color=colour, alpha=0.75, linewidth=1.3, zorder=2)

        ax.scatter(sel["PCo1"], sel["PCo2"], s=42, alpha=0.85, color=colour, edgecolor="white", linewidth=0.6, label=f"{lev} (n={len(sel)})", zorder=3)

    ax.axhline(0, color="#BBBBBB", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#BBBBBB", linewidth=0.8, zorder=0)
    ax.set_xlabel(f"PCo1 ({expl[0]:.1%})", fontsize=11)
    ax.set_ylabel(f"PCo2 ({expl[1]:.1%})", fontsize=11)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
    
    path = os.path.join(OUTPUT_DIR, fname)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✅ PNG Saved: {path}")

rg = results[g_col]
subtitle = f"{metric_label} | PERMANOVA R2={rg['R2']:.3f}, p={rg['p_value']:.4f}"

scatter_png(g_col, f"{args.prefix}_PCoA_by_{g_col}.png", f"{FEATURE_LABEL} Sample Ordination by {g_col} (PCoA)", subtitle)

fig_html = px.scatter(
    ord_df, x="PCo1", y="PCo2", color=g_col,
    hover_name=ord_df.index,
    title=f"{FEATURE_LABEL} Sample Ordination by {g_col} (PCoA)",
    labels={"PCo1": f"PCo1 ({expl[0]:.1%})", "PCo2": f"PCo2 ({expl[1]:.1%})"},
    template="plotly_white"
)
html_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_PCoA_by_{g_col}.html")
fig_html.write_html(html_path)
print(f"✅ HTML Saved: {html_path}")

# =========================================================
# 8. WRITE STATISTICAL SUMMARY REPORTS
# =========================================================
ord_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_pcoa_coordinates.tsv")
ord_df.to_csv(ord_path, sep="\t")

stats_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_permanova.json")
payload = dict(
    matrix=os.path.abspath(matrix_file),
    metadata=os.path.abspath(meta_file),
    metric=args.metric,
    samples=int(df.shape[0]),
    features=int(df.shape[1]),
    permanova=results
)
with open(stats_path, "w") as fh:
    json.dump(payload, fh, indent=2, default=str)

if HAS_MC:
    mc.record(matrix=os.path.abspath(matrix_file), metadata=os.path.abspath(meta_file), output_dir=os.path.abspath(OUTPUT_DIR), prefix=args.prefix)
    mc.write_manifest(OUTPUT_DIR, args.prefix)

print("🎉 Ordination and PERMANOVA execution completed successfully!")
