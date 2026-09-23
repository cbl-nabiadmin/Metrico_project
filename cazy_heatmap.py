#!/usr/bin/env python3
"""
CAZy Class Relative Abundance Heatmap Pipeline
----------------------------------------------
This script processes a Carbohydrate-Active enZYmes (CAZy) class count matrix 
and sample metadata to generate:
1. A static, publication-ready PNG heatmap (via Seaborn Clustermap)
2. An interactive HTML heatmap with a metadata host track (via Plotly)
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

# Configure Matplotlib backend for headless environments (e.g., Docker containers or HPC clusters)
import matplotlib
matplotlib.use('Agg')  # Prevents GUI display errors when saving figures in non-interactive sessions
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# =========================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =========================================================
parser = argparse.ArgumentParser(
    description="CAZy Class Relative Abundance Heatmap Pipeline",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-i", "--input", 
    type=str, 
    default="cazy_class_matrix.tsv",
    help="Filename or relative path of the CAZy class matrix (inside input_file/)"
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
    default="CAZy_Classes_Relative_Abundance",
    help="Prefix name for saved output files"
)
parser.add_argument(
    "-n", "--top-n", 
    type=int, 
    default=30,
    help="Number of top abundant CAZy classes to plot"
)

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

cazy_matrix_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
metadata_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

print(f"📂 Using CAZy Matrix File: {cazy_matrix_file}")
print(f"📂 Using Metadata File:    {metadata_file}")
print(f"📂 Saving Outputs To:      {OUTPUT_DIR}")
print(f"📊 Displaying Top Classes: {args.top_n}")

# =========================================================
# STEP 1: READ MATRIX DATA
# =========================================================
print("Reading CAZy matrix data...")
if not os.path.exists(cazy_matrix_file):
    print(f"❌ CRITICAL ERROR: Input file not found at {cazy_matrix_file}")
    sys.exit(1)

try:
    df_mat = pd.read_csv(cazy_matrix_file, sep=None, engine="python", index_col=0)
except Exception:
    df_mat = pd.read_csv(cazy_matrix_file, sep="\t", index_col=0)

df_mat.index = df_mat.index.astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)

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

meta = meta.set_index("SampleID")[["Host"]]
target_hosts = list(meta["Host"].unique())

# =========================================================
# STEP 3: MATCH MATRIX AND METADATA SAMPLES
# =========================================================
common = list(set(df_mat.index).intersection(set(meta.index)))
print(f"📊 DIAGNOSTIC: Overlapping samples matched: {len(common)}")

if len(common) == 0:
    print("❌ CRITICAL ERROR: Zero sample overlap between CAZy matrix and Metadata!")
    sys.exit(1)

df_mat = df_mat.loc[common]
meta = meta.loc[common]

# =========================================================
# STEP 4: CONVERT VALUES TO NUMERIC
# =========================================================
df_mat = df_mat.apply(pd.to_numeric, errors="coerce").fillna(0)

# =========================================================
# STEP 5: CALCULATE RELATIVE ABUNDANCE (ROW-WISE PERCENTAGE)
# =========================================================
row_totals = df_mat.sum(axis=1)
row_totals[row_totals == 0] = 1.0

df_rel = df_mat.div(row_totals, axis=0) * 100

# =========================================================
# STEP 6: SELECT TOP ABUNDANT CAZy CLASSES
# =========================================================
top_n = min(args.top_n, df_rel.shape[1])
top_cazy = df_rel.sum(axis=0).sort_values(ascending=False).head(top_n).index
df_rel = df_rel[top_cazy]

# =========================================================
# STEP 7: ORDER SAMPLES BY HOST GROUP
# =========================================================
unique_target_hosts = sorted([h for h in target_hosts if h in meta["Host"].unique()])
meta["Host"] = pd.Categorical(meta["Host"], categories=unique_target_hosts, ordered=True)
meta = meta.sort_values("Host")

ordered_samples = meta.index
df_rel = df_rel.loc[ordered_samples]

# =========================================================
# STEP 8: ASSIGN METADATA ANNOTATION COLORS
# =========================================================
host_levels = list(meta["Host"].cat.categories)
palette = sns.color_palette("Set2", n_colors=max(1, len(host_levels)))
host_colors = dict(zip(host_levels, palette))
sample_colors = meta["Host"].astype(str).map(host_colors)

# =========================================================
# STEP 9: DEFINE HEATMAP COLORMAP AND QUANTILE SCALING
# =========================================================
heat_colors = ["#ffffff", "#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"]
cmap = LinearSegmentedColormap.from_list("cazy_heat", heat_colors, N=100)

max_val = float(np.quantile(df_rel.values, 0.95)) if df_rel.size > 0 else 1.0

df_plot = df_rel.T
formatted_cols = [col.replace("_", "\n") for col in df_plot.columns]

# =========================================================
# STEP 10: GENERATE AND SAVE STATIC PNG HEATMAP
# =========================================================
print("Generating static PNG heatmap...")

fig_width = max(12, df_plot.shape[1] * 0.18 + 2.0)
fig_height = max(7, df_plot.shape[0] * 0.30 + 1.5)

# 1. Render heatmap with colorbar positioned at top-left
g = sns.clustermap(
    df_plot,
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

# 2. Safely fetch colorbar axis (supports both old and new Seaborn versions)
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
ax.set_xticklabels(formatted_cols, rotation=90, fontweight="bold", fontsize=6)

# Force Y-axis tick labels to the LEFT side to match Plotly
ax.yaxis.tick_left()
ax.yaxis.set_label_position("left")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontweight="bold", fontsize=8)

ax.set_title(f"Top {top_n} CAZy Classes (Relative Abundance %)", fontsize=13, fontweight="bold", pad=25)

# Place metadata color legend cleanly on the right side
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

png_filename = f"{args.prefix}_CAZy_class.png"
png_output_path = os.path.join(OUTPUT_DIR, png_filename)
plt.savefig(png_output_path, dpi=600, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ PNG Saved: {png_output_path}")

# =========================================================
# STEP 11: GENERATE INTERACTIVE HTML HEATMAP (PLOTLY)
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

if len(host_levels) > 1:
    track_colorscale = [[i / (len(host_levels) - 1), col] for i, col in enumerate(host_hex_palette)]
else:
    track_colorscale = [[0.0, host_hex_palette[0]], [1.0, host_hex_palette[0]]]

fig_html.add_trace(
    go.Heatmap(
        z=[sample_hosts_num],
        x=df_plot.columns,
        y=["Host"],
        colorscale=track_colorscale,
        showscale=False,
        hoverinfo="text",
        text=[[f"Sample: {s}<br>Host: {h}" for s, h in zip(df_plot.columns, sample_hosts_str)]]
    ),
    row=1, col=1
)

hover_text_matrix = []
for cazy_class in df_plot.index:
    row_text = []
    for sample_id in df_plot.columns:
        val = df_plot.loc[cazy_class, sample_id]
        
        if isinstance(val, pd.Series):
            rel_ab = float(val.iloc[0])
        else:
            rel_ab = float(val)

        host_val = meta.loc[sample_id, "Host"]
        if isinstance(host_val, pd.Series):
            host_name = str(host_val.iloc[0])
        else:
            host_name = str(host_val)

        row_text.append(f"SampleID: {sample_id}<br>Host: {host_name}<br>CAZy Class: {cazy_class}<br>Rel Abundance: {rel_ab:.2f}%")
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

fig_html.update_layout(
    title=f"Top {top_n} CAZy Classes (Relative Abundance %) - Interactive View",
    xaxis2=dict(tickangle=90, tickfont=dict(size=8, color="black")),
    # Reverse autorange to keep top-to-bottom row ordering consistent with PNG
    yaxis2=dict(tickfont=dict(size=9, color="black"), autorange="reversed"),
    template="plotly_white",
    height=max(600, df_plot.shape[0] * 20)
)

html_filename = f"{args.prefix}_CAZy_class.html"
html_output_path = os.path.join(OUTPUT_DIR, html_filename)
fig_html.write_html(html_output_path)
print(f"✅ HTML Saved: {html_output_path}")

print("🎉 All tasks completed successfully!")
