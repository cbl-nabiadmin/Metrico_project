#!/usr/bin/env python3
"""
COG Categories Relative Abundance Pipeline
------------------------------------------
This script calculates the sample-level and host-level relative abundance (%) 
for Clusters of Orthologous Groups (COG) functional categories. 

For every COG category, it generates:
1. Static publication-ready PNG bar plots with overlaid individual sample data points (Seaborn)
2. Interactive HTML bar plots with overlaid individual sample data points (Plotly)
3. An all-in-one interactive HTML summary heatmap across all COG categories
4. A static publication-ready 300 DPI PNG summary heatmap across all COG categories
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np

# Configure Matplotlib backend for non-interactive / headless execution (e.g. Docker, HPC pipelines)
import matplotlib
matplotlib.use('Agg')  # Prevents GUI display errors when rendering images without an X11 server
import matplotlib.pyplot as plt

import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go

# =========================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =========================================================
parser = argparse.ArgumentParser(
    description="COG Categories Relative Abundance Pipeline",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument(
    "-i", "--input", 
    type=str, 
    default="all_cog_matrix.csv",
    help="Filename or relative path of the COG matrix (inside input_file/)"
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
    default="/Metrico_project/output/COG_plots",
    help="Path to the output directory where plots will be saved"
)
parser.add_argument(
    "--prefix", 
    type=str, 
    default="COG",
    help="Prefix name for saved output files"
)

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

# Ensure output directory exists before generating figures
os.makedirs(OUTPUT_DIR, exist_ok=True)

cog_matrix_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
metadata_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

print(f"📂 Using COG Matrix File: {cog_matrix_file}")
print(f"📂 Using Metadata File:    {metadata_file}")
print(f"📂 Saving Outputs To:      {OUTPUT_DIR}")

# Lookup dictionary mapping single-letter COG codes to full functional category descriptions
COG_NAMES = {
    "J": "J - Translation, ribosomal structure & biogenesis",
    "A": "A - RNA processing and modification",
    "K": "K - Transcription",
    "L": "L - Replication, recombination and repair",
    "B": "B - Chromatin structure and dynamics",
    "D": "D - Cell cycle control, cell division, chromosome partitioning",
    "Y": "Y - Nuclear structure",
    "V": "V - Defense mechanisms",
    "T": "T - Signal transduction mechanisms",
    "M": "M - Cell wall/membrane/envelope biogenesis",
    "N": "N - Cell motility",
    "Z": "Z - Cytoskeleton",
    "W": "W - Extracellular structures",
    "U": "U - Intracellular trafficking, secretion, and vesicular transport",
    "O": "O - Posttranslational modification, protein turnover, chaperones",
    "C": "C - Energy production and conversion",
    "G": "G - Carbohydrate transport and metabolism",
    "E": "E - Amino acid transport and metabolism",
    "F": "F - Nucleotide transport and metabolism",
    "I": "I - Lipid transport and metabolism",
    "P": "P - Inorganic ion transport and metabolism",
    "Q": "Q - Secondary metabolites biosynthesis, transport and catabolism",
    "R": "R - General function prediction only",
    "S": "S - Function unknown"
}

# =========================================================
# 3. READ AND CLEAN COG MATRIX DATA
# =========================================================
print("Reading COG matrix...")
if not os.path.exists(cog_matrix_file):
    print(f"❌ CRITICAL ERROR: COG file not found at {cog_matrix_file}")
    sys.exit(1)

try:
    mat = pd.read_csv(cog_matrix_file, sep=None, engine='python')
except Exception:
    mat = pd.read_csv(cog_matrix_file, sep="\t")

# Clean column headers
mat.columns = mat.columns.str.strip()

# Ensure the first column is named 'SampleID'
if "SampleID" not in mat.columns:
    mat = mat.rename(columns={mat.columns[0]: "SampleID"})

# Standardize sample names by removing tool suffixes and leading/trailing whitespace
mat["SampleID"] = mat["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)

# Detect all COG functional category columns (all columns except 'SampleID')
cog_cols = [c for c in mat.columns if c.lower() != "sampleid"]
print(f"📊 DIAGNOSTIC: Total COG categories detected ({len(cog_cols)}): {cog_cols}")

# Convert count values to numeric and handle missing values safely
mat[cog_cols] = mat[cog_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

# =========================================================
# 4. READ AND CLEAN METADATA
# =========================================================
print("Reading metadata...")
if not os.path.exists(metadata_file):
    print(f"❌ CRITICAL ERROR: Metadata file not found at {metadata_file}")
    sys.exit(1)

meta = pd.read_csv(metadata_file)
meta.columns = meta.columns.str.strip()

# Dynamically locate Sample and Host columns regardless of exact header casing
s_col = [c for c in meta.columns if "sample" in c.lower()][0] if any("sample" in c.lower() for c in meta.columns) else meta.columns[0]
h_col = [c for c in meta.columns if "host" in c.lower()][0] if any("host" in c.lower() for c in meta.columns) else meta.columns[1]

meta = meta.rename(columns={s_col: "SampleID", h_col: "Host"})
meta["SampleID"] = meta["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta["Host"] = meta["Host"].astype(str).str.strip()

# Exclude accidental header rows
meta = meta[meta["SampleID"].str.lower() != "sampleid"].copy()

# =========================================================
# 5. MERGE DATASETS & CALCULATE TOTAL SAMPLE-LEVEL GENE COUNTS
# =========================================================
merged = pd.merge(mat, meta[["SampleID", "Host"]], on="SampleID", how="inner")
print(f"📊 DIAGNOSTIC: Overlapping samples matched: {len(merged)}")

if merged.empty:
    print("❌ CRITICAL ERROR: Zero overlapping SampleIDs found!")
    sys.exit(1)

# Sum total COG gene occurrences per sample across ALL categories to derive true denominator
merged["Total_COG_Genes"] = merged[cog_cols].sum(axis=1)
merged["Total_COG_Genes"] = merged["Total_COG_Genes"].replace(0, 1.0)  # Prevent division by zero

# Setup distinct categorical color palette for plant host groups
unique_hosts = sorted(merged["Host"].unique())
palette = sns.color_palette("Set2", n_colors=max(1, len(unique_hosts)))
host_colors = dict(zip(unique_hosts, palette))
host_hex_colors = [matplotlib.colors.to_hex(c) for c in palette]
host_hex_map = dict(zip(unique_hosts, host_hex_colors))

# Apply clean whitegrid background style
plt.style.use("seaborn-v0_8-whitegrid") if "seaborn-v0_8-whitegrid" in plt.style.available else plt.style.use("default")

# Container array to hold host summary metrics for building the final summary heatmap
summary_records = []

# =========================================================
# 6. PROCESS EACH COG CATEGORY INDIVIDUALLY
# =========================================================
print(f"\n🚀 Processing plots for {len(cog_cols)} COG categories...")

for cog in cog_cols:
    # 1. Calculate relative abundance (%) for the current COG category per sample
    rel_col_name = f"RelAbun_{cog}"
    merged[rel_col_name] = (merged[cog] / merged["Total_COG_Genes"]) * 100

    # 2. Group by Host to compute mean relative abundance and Standard Error of the Mean (SEM)
    host_df = (
        merged.groupby("Host")[rel_col_name]
        .agg(
            Relative_Abundance="mean",
            SEM="sem"
        )
        .reset_index()
    )

    host_df = host_df.sort_values("Host").reset_index(drop=True)

    # 3. Store mean host values for master summary heatmap
    for _, row in host_df.iterrows():
        summary_records.append({
            "Host": row["Host"],
            "COG": cog,
            "Relative_Abundance": row["Relative_Abundance"]
        })

    # Prepare descriptive title text for plots
    cog_title_desc = COG_NAMES.get(cog, f"Category {cog}")
    title_text = f"COG Category: {cog_title_desc}"

    # Calculate dynamic upper Y-limit based on peak individual sample points
    max_y = merged[rel_col_name].max()
    y_upper = max_y * 1.25 if max_y > 0 else 1.0

    # -----------------------------------------------------
    # A. Generate Static PNG Bar Plot (Matplotlib / Seaborn)
    # -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 6))

    # 1. Draw host mean abundance bars
    sns.barplot(
        data=host_df,
        x="Host",
        y="Relative_Abundance",
        hue="Host",
        palette=host_colors,
        width=0.7,
        ax=ax,
        legend=False,
        alpha=0.65
    )

    # 2. Overlay individual sample data points
    sns.stripplot(
        data=merged,
        x="Host",
        y=rel_col_name,
        hue="Host",
        palette=host_colors,
        size=6.5,
        jitter=0.2,
        edgecolor="black",
        linewidth=1,
        alpha=0.85,
        ax=ax,
        legend=False
    )

    # Axis labels and titles with added padding
    ax.set_title(title_text, fontsize=14, fontweight="bold", pad=20)
    ax.set_xlabel("Host", fontsize=12, fontweight="bold", labelpad=12)
    ax.set_ylabel("Relative Abundance (% of Total COGs)", fontsize=12, fontweight="bold", labelpad=15)

    # Axis tick configuration
    ax.tick_params(axis="x", labelsize=10, pad=5)
    ax.tick_params(axis="y", labelsize=10, pad=5)
    plt.xticks(rotation=45, ha="right", fontweight="bold")
    plt.yticks(fontweight="bold")

    ax.set_ylim(0, y_upper)

    # Clean borders
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0.02, 0.05, 0.98, 0.95])

    # Save PNG figure
    png_filename = f"{args.prefix}_{cog}_relative_abundance.png"
    png_path = os.path.join(OUTPUT_DIR, png_filename)
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close()

    # -----------------------------------------------------
    # B. Generate Interactive HTML Bar Plot (Plotly)
    # -----------------------------------------------------
    # 1. Base bar plot displaying mean host values
    fig_html = px.bar(
        host_df,
        x="Host",
        y="Relative_Abundance",
        color="Host",
        color_discrete_sequence=host_hex_colors,
        title=title_text,
        labels={"Relative_Abundance": "Relative Abundance (% of COGs)", "Host": "Host"},
        text_auto=".2f",
        opacity=0.65
    )

    # 2. Overlay individual interactive sample data points per host
    for h in unique_hosts:
        h_sub = merged[merged["Host"] == h]
        fig_html.add_trace(
            go.Box(
                x=h_sub["Host"],
                y=h_sub[rel_col_name],
                name=h,
                boxpoints='all',
                jitter=0.4,
                pointpos=0,
                fillcolor='rgba(0,0,0,0)',
                line=dict(color='rgba(0,0,0,0)'),
                marker=dict(
                    size=7,
                    color=host_hex_map[h],
                    line=dict(width=1, color='black')
                ),
                hovertext=h_sub["SampleID"],
                hovertemplate="<b>Sample:</b> %{hovertext}<br><b>Relative Abun:</b> %{y:.2f}%<extra></extra>",
                showlegend=False
            )
        )

    fig_html.update_layout(
        xaxis=dict(tickangle=-45, tickfont=dict(size=11, color="black")),
        yaxis=dict(tickfont=dict(size=11, color="black"), range=[0, y_upper]),
        font=dict(size=12),
        template="plotly_white",
        showlegend=False
    )

    # Save HTML output file
    html_filename = f"{args.prefix}_{cog}_relative_abundance.html"
    html_path = os.path.join(OUTPUT_DIR, html_filename)
    fig_html.write_html(html_path)

    print(f"   [✓] COG '{cog}' saved -> PNG & HTML")

# =========================================================
# 7. GENERATE MASTER SUMMARY HEATMAP (ALL COGs ACROSS HOSTS)
# =========================================================
print("\nGenerating All-in-One COG Summary Heatmaps (HTML & PNG)...")
summary_df = pd.DataFrame(summary_records)

# Map single COG letters to full functional category descriptions (A–Z)
summary_df["COG_Description"] = summary_df["COG"].map(
    lambda c: COG_NAMES.get(c, f"Category {c}")
)

# Reshape summary records into a pivot table matrix (Full COG Description vs Host)
summary_pivot = summary_df.pivot(
    index="COG_Description", 
    columns="Host", 
    values="Relative_Abundance"
).fillna(0)

# ---------------------------------------------------------
# A. Interactive HTML Summary Heatmap (Plotly)
# ---------------------------------------------------------
fig_summary_html = px.imshow(
    summary_pivot,
    labels=dict(x="Host", y="COG Category", color="Relative Abundance (%)"),
    x=summary_pivot.columns,
    y=summary_pivot.index,
    color_continuous_scale="YlGnBu",
    title="Overview: True Mean Relative Abundance (% of Total COGs) Across Hosts"
)

fig_summary_html.update_layout(
    xaxis=dict(tickangle=-45, tickfont=dict(size=11, color="black")),
    yaxis=dict(tickfont=dict(size=10, color="black"), autorange="reversed"),
    margin=dict(l=300, r=50, t=80, b=100),  # Generous 300px left margin for long label strings
    height=max(600, len(summary_pivot) * 25),
    template="plotly_white"
)

summary_html_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_ALL_COG_categories_summary_heatmap.html")
fig_summary_html.write_html(summary_html_path)
print(f"✅ Summary Heatmap (HTML) Saved: {summary_html_path}")

# ---------------------------------------------------------
# B. Static PNG Summary Heatmap (Seaborn / Matplotlib)
# ---------------------------------------------------------
fig_height = max(10, len(summary_pivot) * 0.45)
fig_width = max(12, len(summary_pivot.columns) * 1.5 + 6)

fig_png, ax_png = plt.subplots(figsize=(fig_width, fig_height))

sns.heatmap(
    summary_pivot,
    annot=True,
    fmt=".2f",
    cmap="YlGnBu",
    cbar_kws={'label': 'Relative Abundance (% of Total COGs)'},
    ax=ax_png,
    linewidths=0.5,
    linecolor="lightgrey"
)

ax_png.set_title("Overview: True Mean Relative Abundance (% of Total COGs) Across Hosts", fontsize=14, fontweight="bold", pad=20)
ax_png.set_xlabel("Host", fontsize=12, fontweight="bold", labelpad=12)
ax_png.set_ylabel("COG Category", fontsize=12, fontweight="bold", labelpad=15)

plt.xticks(rotation=45, ha="right", fontweight="bold", fontsize=10)
plt.yticks(rotation=0, fontweight="bold", fontsize=10)

plt.tight_layout()

summary_png_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_ALL_COG_categories_summary_heatmap.png")
plt.savefig(summary_png_path, dpi=300, bbox_inches="tight")
plt.close()

print(f"✅ Summary Heatmap (PNG) Saved:  {summary_png_path}")

print(f"\n🎉 All {len(cog_cols)} COG category plots generated successfully in: {OUTPUT_DIR}")
