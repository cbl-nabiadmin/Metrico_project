#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Enables headless rendering in Docker
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# =========================================================
# 1. COMMAND LINE ARGUMENTS
# =========================================================
parser = argparse.ArgumentParser(
    description="EC Numbers Matrix Relative Abundance Heatmap Pipeline",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-i", "--input", 
    type=str, 
    default="all_ec_matrix.tsv",
    help="Filename or relative path of the EC matrix (inside input_file/)"
)
parser.add_argument(
    "-m", "--metadata", 
    type=str, 
    default="Plant-host.csv",
    help="Filename or relative path of the metadata CSV (inside input_file/)"
)
parser.add_argument(
    "-o", "--output-dir", 
    type=str, 
    default="/Metrico_project/output",
    help="Path to the output directory where plots will be saved"
)
parser.add_argument(
    "--prefix", 
    type=str, 
    default="EC_Numbers_Relative_Abundance",
    help="Prefix name for saved output files"
)
parser.add_argument(
    "-n", "--top-n", 
    type=int, 
    default=30,
    help="Number of top abundant EC numbers to plot"
)

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

ec_matrix_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
metadata_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

print(f"📂 Using EC Matrix File: {ec_matrix_file}")
print(f"📂 Using Metadata File:  {metadata_file}")
print(f"📂 Saving Outputs To:     {OUTPUT_DIR}")
print(f"📊 Displaying Top ECs:   {args.top_n}")

# =========================================================
# 🔷 STEP 1: READ EC MATRIX
# =========================================================
print("Reading EC matrix...")
if not os.path.exists(ec_matrix_file):
    print(f"❌ CRITICAL ERROR: EC matrix file not found at {ec_matrix_file}")
    sys.exit(1)

# Auto-detect delimiter (tsv/csv)
try:
    df = pd.read_csv(ec_matrix_file, sep=None, engine="python", index_col=0)
except Exception:
    df = pd.read_csv(ec_matrix_file, sep="\t", index_col=0)

df = df.loc[~df.index.isin(["EC", "Sample"])]
df = df.loc[:, ~df.columns.isin(["Sample", "EC"])]
df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

# =========================================================
# 🔷 STEP 2: READ & CLEAN METADATA
# =========================================================
print("Reading metadata...")
if not os.path.exists(metadata_file):
    print(f"❌ CRITICAL ERROR: Metadata file not found at {metadata_file}")
    sys.exit(1)

meta = pd.read_csv(metadata_file)
meta.columns = meta.columns.str.strip()

# Dynamic sample and host column matching
s_col = [c for c in meta.columns if "sample" in c.lower()][0] if any("sample" in c.lower() for c in meta.columns) else meta.columns[0]
h_col = [c for c in meta.columns if "host" in c.lower()][0] if any("host" in c.lower() for c in meta.columns) else meta.columns[1]

meta = meta.rename(columns={s_col: "SampleID", h_col: "Host"})
meta = meta[meta["SampleID"].astype(str).str.lower() != "sampleid"].copy()

meta["SampleID"] = meta["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta["Host"] = meta["Host"].astype(str).str.strip()
meta = meta.set_index("SampleID")

df.columns = df.columns.astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)

# =========================================================
# 🔷 STEP 3: MATCH DATA + METADATA
# =========================================================
common = list(set(df.columns).intersection(set(meta.index)))
print(f"📊 DIAGNOSTIC: Overlapping samples matched: {len(common)}")

if len(common) == 0:
    print("❌ CRITICAL ERROR: Zero sample overlap between EC matrix and Metadata!")
    sys.exit(1)

df = df[common]
meta = meta.loc[common]

# =========================================================
# 🔷 STEP 4: FILTER LOW ABUNDANCE & SELECT TOP ECs
# =========================================================
df = df[df.sum(axis=1) > 0]
top_n = min(args.top_n, df.shape[0])
top_ec = df.sum(axis=1).sort_values(ascending=False).head(top_n).index
df_top = df.loc[top_ec]

# =========================================================
# 🔷 STEP 5: RELATIVE ABUNDANCE (%)
# =========================================================
col_sums = df_top.sum(axis=0)
col_sums[col_sums == 0] = 1.0
df_rel = df_top.div(col_sums, axis=1).fillna(0).replace([np.inf, -np.inf], 0) * 100

# =========================================================
# 🔷 STEP 6: ORDER SAMPLES BY HOST
# =========================================================
meta["Host"] = meta["Host"].astype("category")
meta = meta.sort_values("Host")

ordered_samples = meta.index
df_rel = df_rel[ordered_samples]

# =========================================================
# 🔷 STEP 7: COLORS
# =========================================================
heat_colors = ["#ffffff", "#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"]
cmap = LinearSegmentedColormap.from_list("custom_heat", heat_colors, N=100)

host_levels = list(meta["Host"].cat.categories)
palette = sns.color_palette("Set2", n_colors=max(1, len(host_levels)))
host_colors = dict(zip(host_levels, palette))

# Convert host series to string before mapping
sample_colors = meta["Host"].astype(str).map(host_colors)

max_val = float(np.quantile(df_rel.values, 0.95)) if df_rel.size > 0 else 1.0
formatted_cols = [col.replace("_", "\n") for col in df_rel.columns]

# =========================================================
# 🔷 STEP 8: SAVE STATIC PNG PLOT
# =========================================================
print("Generating static PNG heatmap...")

fig_width = max(12, df_rel.shape[1] * 0.18 + 2.0)
fig_height = max(7, df_rel.shape[0] * 0.30 + 1.5)

# 1. Render heatmap with explicit top-left colorbar positioning
g = sns.clustermap(
    df_rel,
    cmap=cmap,
    vmax=max_val if max_val > 0 else 1.0,
    vmin=0,
    col_colors=sample_colors,
    row_cluster=False,
    col_cluster=False,
    figsize=(fig_width, fig_height),
    cbar_pos=(0.02, 0.80, 0.03, 0.12),  # Position [left, bottom, width, height]
    cbar_kws={"label": ""},              # Clear default vertical label string
    linewidths=0,
    yticklabels=True,
    xticklabels=formatted_cols
)

# 2. Safely fetch colorbar axis across Seaborn versions
cbar_axis = getattr(g, "ax_cbar", getattr(g, "cbar_ax", None))

# 3. Add horizontal title directly ON TOP of the color scale legend
if cbar_axis is not None:
    cbar_axis.set_title(
        "Relative Abundance (%)", 
        fontsize=9, 
        fontweight="bold", 
        pad=10, 
        ha="center"  # Centered horizontally over the color scale bar
    )
    cbar_axis.tick_params(labelsize=8)

ax = g.ax_heatmap
ax.set_xticklabels(formatted_cols, rotation=90, fontweight="bold", fontsize=7)

# Force y-axis labels to the LEFT side in Seaborn to match Plotly
ax.yaxis.tick_left()
ax.yaxis.set_label_position("left")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontweight="bold", fontsize=7)

ax.set_title(f"Top {top_n} EC Numbers Relative Abundance (EggNOG)", fontsize=14, fontweight="bold", pad=25)

legend_patches = [mpatches.Patch(color=color, label=host) for host, color in host_colors.items()]
if g.ax_col_dendrogram:
    g.ax_col_dendrogram.legend(
        handles=legend_patches,
        title="Host",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=9,
        title_fontsize=10
    )

png_filename = f"{args.prefix}_EC_matrix.png"
png_output_path = os.path.join(OUTPUT_DIR, png_filename)
plt.savefig(
    png_output_path,
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)
plt.close()
print(f"✅ PNG Saved: {png_output_path}")

# =========================================================
# 🔷 STEP 9: SAVE INTERACTIVE HTML HEATMAP WITH HOST TRACK
# =========================================================
print("Generating interactive HTML heatmap...")

host_to_int = {h: idx for idx, h in enumerate(host_levels)}
sample_hosts_str = [str(meta.loc[s, "Host"]) for s in df_rel.columns]
sample_hosts_num = [host_to_int[h] for h in sample_hosts_str]

fig_html = make_subplots(
    rows=2, cols=1, 
    row_heights=[0.05, 0.95], 
    shared_xaxes=True, 
    vertical_spacing=0.02
)

host_hex_palette = [sns.color_palette("Set2", n_colors=max(1, len(host_levels))).as_hex()[i] for i in range(len(host_levels))]

# 1. Top Track: Host Color Metadata
fig_html.add_trace(
    go.Heatmap(
        z=[sample_hosts_num],
        x=df_rel.columns,
        y=["Host"],
        colorscale=[[i / max(1, len(host_levels) - 1), col] for i, col in enumerate(host_hex_palette)] if len(host_levels) > 1 else [[0, host_hex_palette[0]], [1, host_hex_palette[0]]],
        showscale=False,
        hoverinfo="text",
        text=[[f"Sample: {s}<br>Host: {h}" for s, h in zip(df_rel.columns, sample_hosts_str)]]
    ),
    row=1, col=1
)

# 2. Main Heatmap: EC Relative Abundance
hover_text_matrix = []
for ec_num in df_rel.index:
    row_text = []
    for sample_id in df_rel.columns:
        rel_ab = df_rel.loc[ec_num, sample_id]
        host_name = meta.loc[sample_id, "Host"]
        row_text.append(f"SampleID: {sample_id}<br>Host: {host_name}<br>EC Number: {ec_num}<br>Rel Abundance: {rel_ab:.2f}%")
    hover_text_matrix.append(row_text)

fig_html.add_trace(
    go.Heatmap(
        z=df_rel.values,
        x=df_rel.columns,
        y=df_rel.index,
        colorscale=heat_colors,
        zmin=0,
        zmax=max_val if max_val > 0 else 1.0,
        colorbar=dict(title="Relative Abundance (%)"),
        hoverinfo="text",
        text=hover_text_matrix
    ),
    row=2, col=1
)

fig_html.update_layout(
    title=f"Top {top_n} EC Numbers Relative Abundance (%) - Interactive View",
    xaxis2=dict(tickangle=90, tickfont=dict(size=8, color="black")),
    # Reverse autorange so row 0 (most abundant EC) starts at the top
    yaxis2=dict(tickfont=dict(size=8, color="black"), autorange="reversed"),
    template="plotly_white",
    height=max(600, df_rel.shape[0] * 20)
)

html_filename = f"{args.prefix}_EC_matrix.html"
html_output_path = os.path.join(OUTPUT_DIR, html_filename)
fig_html.write_html(html_output_path)
print(f"✅ HTML Saved: {html_output_path}")

print("🎉 All tasks completed successfully!")
