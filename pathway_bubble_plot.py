#!/usr/bin/env python3
"""
Peak Taxon Pathway Driver Bubble Plot Pipeline
------------------------------------------------------------------------
Processes HUMAnN pathway abundance data and metadata. Drops UNMAPPED and 
UNINTEGRATED entries completely, preserves unclassified taxon reads for 
specific pathways, identifies the peak contributing taxon driver per 
(Sample, Pathway) pair, and plots top pathways across distinct plant hosts.
"""

import os
import re
import sys
import argparse
import numpy as np
import pandas as pd

# Configure Matplotlib backend for headless rendering
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# =========================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =========================================================
parser = argparse.ArgumentParser(
    description="Peak Taxon Pathway Driver Bubble Plot Pipeline",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("-i", "--input", type=str, default="all-pathabundance_cpm.tsv", help="Pathway abundance file")
parser.add_argument("-m", "--metadata", type=str, default="Plant-host.csv", help="Metadata CSV")
parser.add_argument("--pathway-meta", type=str, default="MetaCyc-pathways.txt", help="Pathway ontology file")
parser.add_argument("-o", "--output-dir", type=str, default="/Metrico_project/output", help="Output directory")
parser.add_argument("--prefix", type=str, default="pathway_bubble_plot", help="Prefix for output files")

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

abund_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
meta_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

print(f"📂 Using Abundance File:    {abund_file}")
print(f"📂 Using Metadata File:     {meta_file}")
print(f"📂 Saving Outputs To:        {OUTPUT_DIR}")

# =========================================================
# 3. READ METADATA & PATHWAY ABUNDANCE MATRIX DATA
# =========================================================
if not os.path.exists(abund_file) or not os.path.exists(meta_file):
    print("❌ CRITICAL ERROR: Input files not found!")
    sys.exit(1)

meta = pd.read_csv(meta_file)
meta.columns = meta.columns.str.strip()

sample_col = [c for c in meta.columns if "sample" in c.lower()]
meta_sample_col = sample_col[0] if sample_col else meta.columns[0]

host_col = [c for c in meta.columns if "host" in c.lower()]
meta_host_col = host_col[0] if host_col else meta.columns[1]

meta = meta.rename(columns={meta_sample_col: "SampleID", meta_host_col: "Host"})

meta["SampleID"] = meta["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta["Host"] = meta["Host"].astype(str).str.strip()

df_abund = pd.read_csv(abund_file, sep="\t")
df_abund.rename(columns={df_abund.columns[0]: "Pathway"}, inplace=True)
df_abund["Pathway"] = df_abund["Pathway"].astype(str).str.strip()

raw_headers = [str(c).strip() for c in df_abund.columns if c != "Pathway"]

def normalize_id(s):
    return re.sub(r"[^a-zA-Z0-9]", "", str(s)).lower()

meta_id_map = {normalize_id(m_id): m_id for m_id in meta["SampleID"].dropna().unique()}

id_mapping = {}
for h in raw_headers:
    h_clean = re.sub(r"(_merged)?_Abundance.*$", "", h, flags=re.IGNORECASE)
    h_clean = re.sub(r"_CPM.*$", "", h_clean, flags=re.IGNORECASE)
    h_clean = re.sub(r"_eggnog.*$", "", h_clean, flags=re.IGNORECASE)
    
    h_norm = normalize_id(h_clean)
    
    if h_norm in meta_id_map:
        id_mapping[h] = meta_id_map[h_norm]
    else:
        matched = None
        for norm_meta_id, orig_meta_id in meta_id_map.items():
            if norm_meta_id in h_norm or h_norm in norm_meta_id:
                matched = orig_meta_id
                break
        id_mapping[h] = matched if matched else h_clean

print(f"📊 DIAGNOSTIC: Header Mapping Result -> {id_mapping}")

# =========================================================
# 4. STRICTLY DROP UNMAPPED/UNINTEGRATED & PROCESS STRATIFIED TAXA
# =========================================================
print("Processing stratified taxon drivers...")

df_long = df_abund.melt(id_vars=["Pathway"], var_name="RawHeader", value_name="CPM")
df_long["SampleID"] = df_long["RawHeader"].astype(str).map(id_mapping)
df_long["CPM"] = pd.to_numeric(df_long["CPM"], errors="coerce").fillna(0)
df_long = df_long[df_long["CPM"] > 0]

df_long = df_long[
    ~df_long["Pathway"].str.upper().str.startswith("UNMAPPED") &
    ~df_long["Pathway"].str.upper().str.startswith("UNINTEGRATED") &
    ~df_long["Pathway"].str.upper().str.startswith("UNGROUPED")
].copy()

is_stratified = df_long["Pathway"].str.contains(r"\|", regex=True)

if is_stratified.any():
    df_taxa = df_long[is_stratified].copy()
    
    split_df = df_taxa["Pathway"].str.split("|", n=1, expand=True)
    df_taxa["Pathway"] = split_df[0]
    df_taxa["Taxon"] = split_df[1].fillna("unclassified")
    
    df_taxa = df_taxa[~df_taxa["Taxon"].str.lower().str.contains("unclassified", na=False)].copy()
    
    df_taxa["Taxon"] = df_taxa["Taxon"].str.replace(r"^.*\.s__", "", regex=True)
    df_taxa["Taxon"] = df_taxa["Taxon"].str.replace(r"^.*\.g__", "", regex=True)
    df_taxa["Taxon"] = df_taxa["Taxon"].str.replace(r"^g__", "", regex=True)
    df_taxa["Taxon"] = df_taxa["Taxon"].str.replace(r"^s__", "", regex=True)
    df_taxa["Taxon"] = df_taxa["Taxon"].str.replace("_", " ")
    
    if not df_taxa.empty:
        df_long = df_taxa
    else:
        print("⚠️ Warning: No classified species found in stratified lines. Falling back to community total.")
        df_long["Taxon"] = "Community Total"
else:
    df_long["Taxon"] = "Community Total"

df_merged = pd.merge(df_long, meta, on="SampleID", how="inner")
df_merged["Pathway_Clean"] = df_merged["Pathway"].str.replace(r"^.*?: ", "", regex=True)

print(f"📊 DIAGNOSTIC: Total non-zero functional rows matched: {len(df_merged)}")
print(f"📊 DIAGNOSTIC: Unique Taxa Found: {df_merged['Taxon'].nunique()} -> {df_merged['Taxon'].unique()[:5].tolist()}")

if df_merged.empty:
    print("❌ CRITICAL ERROR: Zero matching rows after filtering!")
    sys.exit(1)

top_pathways = df_merged.groupby("Pathway_Clean")["CPM"].sum().nlargest(25).index.tolist()
df_filtered = df_merged[df_merged["Pathway_Clean"].isin(top_pathways)].copy()

# =========================================================
# 5. ISOLATE EXACTLY ONE PEAK CONTRIBUTING TAXON PER PAIR
# =========================================================
print("Extracting top contributing taxon driver for each (Sample, Pathway) pair...")

df_peak = (
    df_filtered.sort_values(by="CPM", ascending=False)
    .groupby(["SampleID", "Pathway_Clean"])
    .first()
    .reset_index()
)

df_peak = df_peak.sort_values(["Host", "SampleID"])
ordered_samples = df_peak["SampleID"].unique()

df_peak["Pathway_Label"] = df_peak["Pathway_Clean"].apply(lambda x: x[:40] + "..." if len(x) > 42 else x)
df_peak["Taxon"] = df_peak["Taxon"].fillna("Unclassified Driver").astype(str)

hosts = sorted([str(h) for h in df_peak["Host"].unique() if pd.notna(h)])

# Ensure colors are converted to HEX so both Seaborn & Plotly use identical colors
host_palette_raw = sns.color_palette("Set2", n_colors=max(1, len(hosts)))
host_hex_colors = [mcolors.to_hex(c) for c in host_palette_raw]
host_colors = dict(zip(hosts, host_hex_colors))

top_taxa = sorted([str(t) for t in df_peak["Taxon"].unique() if pd.notna(t)])
taxon_palette_raw = sns.color_palette("husl", n_colors=max(1, len(top_taxa)))
taxon_hex_colors = [mcolors.to_hex(c) for c in taxon_palette_raw]
taxon_colors = dict(zip(top_taxa, taxon_hex_colors))

# =========================================================
# 6. GENERATE STATIC PNG BUBBLE PLOT
# =========================================================
print("Generating static PNG plot...")

num_pathways = len(df_peak["Pathway_Label"].unique())
num_samples = len(ordered_samples)

# DYNAMIC CANVAS RESIZING: Prevents excessive whitespace
fig_width = max(14, num_samples * 0.45 + 6)
fig_height = max(8, num_pathways * 0.35 + 4)

fig = plt.figure(figsize=(fig_width, fig_height))

gs_main = gridspec.GridSpec(2, 1, height_ratios=[max(4, num_pathways * 0.3), 2.2], hspace=0.25)
gs_top = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_main[0], height_ratios=[0.04, 0.96], hspace=0.02)

ax_host = fig.add_subplot(gs_top[0])
ax_main = fig.add_subplot(gs_top[1])

# Host Track
sample_x_pos = {samp: i for i, samp in enumerate(ordered_samples)}
df_peak["x_pos"] = df_peak["SampleID"].map(sample_x_pos)

sample_host_map = df_peak.groupby("SampleID")["Host"].first()
for samp, x_idx in sample_x_pos.items():
    h = sample_host_map.get(samp)
    ax_host.add_patch(mpatches.Rectangle((x_idx - 0.5, 0), 1, 1, color=host_colors.get(h, "#CCCCCC")))

ax_host.set_xlim(-0.8, len(ordered_samples) - 0.2)
ax_host.set_ylim(0, 1)

ax_host.set_yticks([0.5])
ax_host.set_yticklabels(["Host"], fontweight="bold", fontsize=10)
ax_host.set_xticks([])
ax_host.spines['top'].set_visible(False)
ax_host.spines['right'].set_visible(False)
ax_host.spines['bottom'].set_visible(False)
ax_host.spines['left'].set_visible(False)

# DYNAMIC BUBBLE SCALING FOR PNG
max_cpm = df_peak["CPM"].max() if not df_peak.empty else 1.0
cpm_sizes = [10, 50, 100, 250, 500, 1000]
cpm_sizes = [s for s in cpm_sizes if s <= max_cpm * 1.5]

# Scale marker size dynamically to fit canvas proportions cleanly
marker_area_scale = 350.0 / (max_cpm if max_cpm > 0 else 1.0)
sizes = df_peak["CPM"] * marker_area_scale
colors = [taxon_colors.get(t, "#888888") for t in df_peak["Taxon"]]

ax_main.scatter(
    df_peak["x_pos"], 
    df_peak["Pathway_Label"], 
    s=sizes, 
    c=colors, 
    alpha=0.85, 
    edgecolors="black", 
    linewidths=0.5
)

current_host = None
for i, samp in enumerate(ordered_samples):
    h = sample_host_map.get(samp)
    if current_host is None:
        current_host = h
    elif h != current_host:
        ax_main.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
        current_host = h

ax_main.set_xticks(range(len(ordered_samples)))
ax_main.set_xticklabels(ordered_samples, rotation=90, fontsize=8, fontweight="bold")
ax_main.set_xlim(-0.8, len(ordered_samples) - 0.2)
ax_main.grid(True, linestyle=":", alpha=0.4, axis="y")

ax_main.set_ylabel("MetaCyc Functional Pathways", fontweight="bold", fontsize=11)
ax_main.set_xlabel("Samples (Grouped by Host)", fontweight="bold", fontsize=11)

# LEGEND GRID (Multi-column layout)
gs_bottom = gridspec.GridSpecFromSubplotSpec(
    1, 4, 
    subplot_spec=gs_main[1], 
    width_ratios=[0.18, 0.27, 0.27, 0.28], 
    wspace=0.15
)

ax_leg0 = fig.add_subplot(gs_bottom[0, 0])
ax_leg1 = fig.add_subplot(gs_bottom[0, 1])
ax_leg2 = fig.add_subplot(gs_bottom[0, 2])
ax_leg3 = fig.add_subplot(gs_bottom[0, 3])

for a in [ax_leg0, ax_leg1, ax_leg2, ax_leg3]:
    a.axis("off")

cpm_handles = [
    Line2D([0], [0], marker='o', color='w', label=f"{s} CPM", markerfacecolor='gray', markersize=np.sqrt(s * marker_area_scale)) 
    for s in cpm_sizes
]
leg_cpm = ax_leg0.legend(handles=cpm_handles, title=r"$\bf{Peak\ Abundance\ (CPM)}$", loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=8.5, frameon=False)
ax_leg0.add_artist(leg_cpm)

host_handles = [mpatches.Patch(color=host_colors[h], label=h) for h in hosts]
ax_leg0.legend(handles=host_handles, title=r"$\bf{Plant\ Host}$", loc="upper left", bbox_to_anchor=(0, 0.45), fontsize=8.5, frameon=False)

taxa_splits = np.array_split(top_taxa, 3) if len(top_taxa) >= 3 else [top_taxa]
leg_axes = [ax_leg1, ax_leg2, ax_leg3]
for idx, split_taxa in enumerate(taxa_splits):
    if idx < len(leg_axes):
        handles = [Line2D([0], [0], marker='o', color='w', label=t, markerfacecolor=taxon_colors[t], markersize=6) for t in split_taxa]
        title_text = r"$\bf{Peak\ Taxon\ Driver}$" if idx == 0 else ""
        leg_axes[idx].legend(handles=handles, title=title_text, loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=8.0, frameon=False, ncol=1)

png_output_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_Pathway_bubble.png")
plt.savefig(png_output_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ PNG Saved: {png_output_path}")

# =========================================================
# 7. GENERATE INTERACTIVE HTML PLOT (PLOTLY)
# =========================================================
print("Generating interactive HTML plot...")

host_to_int = {h: idx for idx, h in enumerate(hosts)}
sample_hosts_str = [str(sample_host_map.get(s)) for s in ordered_samples]
sample_hosts_num = [host_to_int[h] for h in sample_hosts_str]

fig_html = make_subplots(
    rows=2, cols=1, 
    row_heights=[0.04, 0.96], 
    shared_xaxes=True, 
    vertical_spacing=0.02
)

if len(hosts) > 1:
    plotly_colorscale = [[i / (len(hosts) - 1), host_colors[hosts[i]]] for i in range(len(hosts))]
else:
    plotly_colorscale = [[0.0, host_colors[hosts[0]]], [1.0, host_colors[hosts[0]]]]

# 1. Top Track: Host Color Metadata Bar
fig_html.add_trace(
    go.Heatmap(
        z=[sample_hosts_num],
        x=ordered_samples,
        y=["Host"],
        colorscale=plotly_colorscale,
        showscale=False,
        hoverinfo="text",
        text=[[f"Sample: {s}<br>Host: {h}" for s, h in zip(ordered_samples, sample_hosts_str)]]
    ),
    row=1, col=1
)

# 2. Main Bubble Plot: Using HEX Colors to Match PNG Exactly
plotly_size_scale = 22.0 / np.sqrt(max_cpm if max_cpm > 0 else 1.0)

for taxon_name in top_taxa:
    df_sub = df_peak[df_peak["Taxon"] == taxon_name]
    if not df_sub.empty:
        fig_html.add_trace(
            go.Scatter(
                x=df_sub["SampleID"],
                y=df_sub["Pathway_Label"],
                mode="markers",
                name=taxon_name,
                marker=dict(
                    size=np.sqrt(df_sub["CPM"]) * plotly_size_scale,
                    color=taxon_colors[taxon_name],  # Uses valid HEX string
                    line=dict(width=0.5, color="black")
                ),
                text=[f"SampleID: {r['SampleID']}<br>Host: {r['Host']}<br>Pathway: {r['Pathway_Clean']}<br>Peak Taxon: {r['Taxon']}<br>Abundance: {r['CPM']:.2f} CPM" for _, r in df_sub.iterrows()],
                hoverinfo="text"
            ),
            row=2, col=1
        )

fig_html.update_layout(
    title="Peak Taxon Driver MetaCyc Pathway Bubble Plot",
    xaxis2=dict(tickangle=90, title="Samples (Grouped by Host)"),
    yaxis2=dict(title="MetaCyc Functional Pathways"),
    template="plotly_white",
    height=max(650, num_pathways * 25 + 200),
    margin=dict(l=300, b=150)
)

html_output_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_Pathway_bubble.html")
fig_html.write_html(html_output_path)
print(f"✅ HTML Saved: {html_output_path}")

print("🎉 Completed successfully!")
