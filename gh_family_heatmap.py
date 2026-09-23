#!/usr/bin/env python3
"""
GH Family Relative Abundance Heatmap Pipeline
----------------------------------------------
This script processes Glycoside Hydrolase (GH) family count data alongside 
sample metadata to compute sample-level relative abundances (%) and generate:
1. A static PNG heatmap (using Seaborn Clustermap) with top host metadata color tracks
2. An interactive HTML heatmap (using Plotly) with custom hover readouts and metadata tracks
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

# Configure Matplotlib for headless rendering (e.g., inside Docker or HPC clusters)
import matplotlib
matplotlib.use('Agg')  # Prevents display errors when running without a graphical interface (X11 server)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# =========================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =========================================================
parser = argparse.ArgumentParser(
    description="GH Family Relative Abundance Heatmap Pipeline",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-i", "--input", 
    type=str, 
    default="all_cazy_count.csv",
    help="Filename or relative path of the CAZy count file (inside input_file/)"
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
    default="GH_Family_Relative_Abundance",
    help="Prefix name for saved output files"
)
parser.add_argument(
    "-n", "--top-n", 
    type=int, 
    default=40,
    help="Number of top abundant GH families to plot"
)

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

cazy_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
metadata_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

print(f"📂 Using CAZy Count File: {cazy_file}")
print(f"📂 Using Metadata File:   {metadata_file}")
print(f"📂 Saving Outputs To:     {OUTPUT_DIR}")
print(f"📊 Displaying Top GHs:    {args.top_n}")

# =========================================================
# STEP 1: READ & PIVOT CAZY COUNT DATA
# =========================================================
print("Reading CAZy counts...")
if not os.path.exists(cazy_file):
    print(f"❌ CRITICAL ERROR: CAZy input file not found at {cazy_file}")
    sys.exit(1)

cazy = pd.read_csv(cazy_file, sep=None, engine='python', header=0)
cazy.columns = cazy.columns.str.strip()

cazy = cazy.rename(columns={
    cazy.columns[0]: "Count",
    cazy.columns[1]: "SampleID",
    cazy.columns[2]: "GH_Family"
})

cazy["SampleID"] = cazy["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
cazy["GH_Family"] = cazy["GH_Family"].astype(str).str.strip()
cazy["Count"] = pd.to_numeric(cazy["Count"], errors="coerce").fillna(0)

cazy = cazy[cazy["SampleID"].str.lower() != "sampleid"].copy()

print(f"📊 DIAGNOSTIC: Total entries: {len(cazy)}")
print(f"📊 DIAGNOSTIC: Unique GH families: {cazy['GH_Family'].nunique()}")

df_mat = cazy.pivot_table(index="SampleID", columns="GH_Family", values="Count", aggfunc="sum", fill_value=0)
print(f"📊 DIAGNOSTIC: Matrix created. Shape: {df_mat.shape} (Samples x GH_Families)")

# =========================================================
# STEP 2: READ AND CLEAN METADATA
# =========================================================
print("Reading metadata...")
if not os.path.exists(metadata_file):
    print(f"❌ CRITICAL ERROR: Metadata file not found at {metadata_file}")
    sys.exit(1)

meta = pd.read_csv(metadata_file)
meta.columns = meta.columns.str.strip()

s_col = [c for c in meta.columns if "sample" in c.lower()][0] if any("sample" in c.lower() for c in meta.columns) else meta.columns[0]
h_col = [c for c in meta.columns if "host" in c.lower()][0] if any("host" in c.lower() for c in meta.columns) else meta.columns[1]

meta = meta.rename(columns={s_col: "SampleID", h_col: "Host"})
meta = meta[meta["SampleID"].astype(str).str.lower() != "sampleid"].copy()
meta["SampleID"] = meta["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta["Host"] = meta["Host"].astype(str).str.strip()

meta = meta.set_index("SampleID")

print(f"📊 DIAGNOSTIC: Total metadata samples: {len(meta)}")

# =========================================================
# STEP 3: SAMPLE INTERSECTION AND MATCHING
# =========================================================
common = list(set(df_mat.index).intersection(set(meta.index)))

print(f"📊 DIAGNOSTIC: Overlapping samples matched: {len(common)}")

if len(common) == 0:
    print("❌ CRITICAL ERROR: Zero overlapping SampleIDs found!")
    sys.exit(1)

df_mat = df_mat.loc[common]
meta = meta.loc[common]
df_mat = df_mat.apply(pd.to_numeric, errors="coerce").fillna(0)

# =========================================================
# STEP 4: RELATIVE ABUNDANCE (%) AND TOP FAMILY SELECTION
# =========================================================
row_sums = df_mat.sum(axis=1)
row_sums[row_sums == 0] = 1.0

df_rel = df_mat.div(row_sums, axis=0).fillna(0) * 100

top_n = min(args.top_n, max(1, df_rel.shape[1]))
top_gh = df_rel.sum(axis=0).sort_values(ascending=False).head(top_n).index
df_rel = df_rel[top_gh]

unique_hosts = sorted(list(meta["Host"].unique()))
meta["Host"] = pd.Categorical(meta["Host"], categories=unique_hosts, ordered=True)
meta = meta.sort_values("Host")

ordered_samples = meta.index
df_rel = df_rel.loc[ordered_samples]

df_plot = df_rel.T

print(f"📊 DIAGNOSTIC: Final Matrix shape for plotting (GH Families x Samples): {df_plot.shape}")

# =========================================================
# STEP 5: COLOR PALETTE & RENDER SETUP
# =========================================================
host_levels = list(meta["Host"].cat.categories)
palette = sns.color_palette("Set2", n_colors=max(1, len(host_levels)))
host_colors = dict(zip(host_levels, palette))

col_colors_df = pd.DataFrame(
    {"Host": meta.loc[df_plot.columns, "Host"].astype(str).map(host_colors)},
    index=df_plot.columns
)

heat_colors = ["#ffffff", "#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"]
cmap = LinearSegmentedColormap.from_list("gh_heat", heat_colors, N=100)

max_val = float(np.quantile(df_plot.values, 0.95)) if df_plot.size > 0 else 1.0

# =========================================================
# STEP 6: DRAW STATIC PNG HEATMAP (SEABORN CLUSTERMAP)
# =========================================================
print("Generating static PNG heatmap...")

fig_width = max(12, df_plot.shape[1] * 0.2)
fig_height = max(8, df_plot.shape[0] * 0.25)

# Render static heatmap without dendrogram clustering to maintain abundance order
g = sns.clustermap(
    df_plot, 
    cmap=cmap,
    vmax=max_val if max_val > 0 else 1.0,
    vmin=0,
    col_colors=col_colors_df, 
    row_cluster=False,                 # Disable row clustering to remove dendrogram
    col_cluster=False,                 # Disable column clustering to preserve host grouping
    figsize=(fig_width, fig_height),
    cbar_pos=(0.02, 0.80, 0.03, 0.12),  # Position top-left colorbar
    cbar_kws={"label": ""},             # Clear vertical label
    linewidths=0, 
    yticklabels=True, 
    xticklabels=True
)

# Safely fetch colorbar axis across Seaborn versions
cbar_axis = getattr(g, "ax_cbar", getattr(g, "cbar_ax", None))

# Add horizontal title directly ON TOP of the color scale box
if cbar_axis is not None:
    cbar_axis.set_title(
        "Relative Abundance (%)", 
        fontsize=9, 
        fontweight="bold", 
        pad=10, 
        ha="center"
    )
    cbar_axis.tick_params(labelsize=8)

ax = g.ax_heatmap
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontweight="bold", fontsize=9)

# Force Y-axis tick labels to the LEFT side
ax.yaxis.tick_left()
ax.yaxis.set_label_position("left")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontweight="bold", fontsize=8)

ax.set_title(f"Top {top_n} GH Family Member Relative Abundance (%)", fontsize=12, fontweight="bold", pad=25)

# Construct and place host metadata color legend
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

png_filename = f"{args.prefix}_GH_Family.png"
png_output_path = os.path.join(OUTPUT_DIR, png_filename)
plt.savefig(png_output_path, dpi=600, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ PNG Saved: {png_output_path}")

# =========================================================
# STEP 7: DRAW INTERACTIVE HTML HEATMAP (PLOTLY)
# =========================================================
print("Generating interactive HTML heatmap...")

host_to_int = {h: idx for idx, h in enumerate(host_levels)}
sample_hosts_str = [str(meta.loc[s, "Host"]) for s in df_plot.columns]
sample_hosts_num = [host_to_int[h] for h in sample_hosts_str]

fig_html = make_subplots(
    rows=2, cols=1, 
    row_heights=[0.05, 0.95], 
    shared_xaxes=True, 
    vertical_spacing=0.02
)

host_hex_palette = [sns.color_palette("Set2", n_colors=max(1, len(host_levels))).as_hex()[i] for i in range(len(host_levels))]

# 1. Top Track: Host Color Metadata Bar
fig_html.add_trace(
    go.Heatmap(
        z=[sample_hosts_num],
        x=df_plot.columns,
        y=["Host"],
        colorscale=[[i / max(1, len(host_levels) - 1), col] for i, col in enumerate(host_hex_palette)] if len(host_levels) > 1 else [[0, host_hex_palette[0]], [1, host_hex_palette[0]]],
        showscale=False,
        hoverinfo="text",
        text=[[f"Sample: {s}<br>Host: {h}" for s, h in zip(df_plot.columns, sample_hosts_str)]]
    ),
    row=1, col=1
)

# 2. Main Heatmap: Relative Abundance (%) Matrix
hover_text_matrix = []
for gh_family in df_plot.index:
    row_text = []
    for sample_id in df_plot.columns:
        rel_ab = df_plot.loc[gh_family, sample_id]
        host_name = meta.loc[sample_id, "Host"]
        row_text.append(f"SampleID: {sample_id}<br>Host: {host_name}<br>GH Family: {gh_family}<br>Relative Abundance: {rel_ab:.2f}%")
    hover_text_matrix.append(row_text)

fig_html.add_trace(
    go.Heatmap(
        z=df_plot.values,
        x=df_plot.columns,
        y=df_plot.index,
        colorscale=heat_colors,
        zmin=0,
        zmax=max_val if max_val > 0 else 1.0,
        colorbar=dict(title="Relative Abundance (%)"),
        hoverinfo="text",
        text=hover_text_matrix
    ),
    row=2, col=1
)

# Configure responsive layout and fonts
fig_html.update_layout(
    title=f"Top {top_n} GH Family Relative Abundance (%) - Interactive View",
    xaxis2=dict(tickangle=90, tickfont=dict(size=8, color="black")),
    # Reverse autorange to keep top-to-bottom row ordering consistent with PNG
    yaxis2=dict(tickfont=dict(size=8, color="black"), autorange="reversed"),
    template="plotly_white",
    height=max(600, df_plot.shape[0] * 18)
)

html_filename = f"{args.prefix}_GH_Family.html"
html_output_path = os.path.join(OUTPUT_DIR, html_filename)
fig_html.write_html(html_output_path)
print(f"✅ HTML Saved: {html_output_path}")

print("🎉 All tasks completed successfully!")
