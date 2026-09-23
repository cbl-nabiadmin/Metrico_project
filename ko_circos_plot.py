#!/usr/bin/env python3
"""
KO Circos Diagram and Module Sharing Pipeline (Proportional Grid Layout)
------------------------------------------------------------------------
Processes long-format KEGG Orthology (KO) IDs and KEGG Module data alongside 
sample metadata. Uses a 2-row multi-column Matplotlib GridSpec layout to 
physically isolate the Circos plot on the top axis and cleanly space all 
legend panels across the bottom grid.
"""

import os
import re
import sys
import argparse
import requests
import numpy as np
import pandas as pd

# Configure Matplotlib backend for headless rendering in server/container environments
import matplotlib
matplotlib.use('Agg')  # Prevents display errors when saving figures without an active X11 display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

import seaborn as sns
from pycirclize import Circos
import plotly.express as px

# =========================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =========================================================
parser = argparse.ArgumentParser(
    description="KO Circos Diagram Pipeline (Proportional Grid Layout)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "-i", "--input-ko", 
    type=str, 
    default="all_sample_ko.tsv",
    help="Filename or relative path of the KO long-format file (inside input_file/)"
)
parser.add_argument(
    "-mod", "--input-module", 
    type=str, 
    default="all_module.tsv",
    help="Filename or relative path of the Module long-format file (inside input_file/)"
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
    default="ko_circos",
    help="Prefix name for saved output files"
)

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

ko_file = args.input_ko if os.path.isabs(args.input_ko) else os.path.join(INPUT_DIR, args.input_ko)
mod_file = args.input_module if os.path.isabs(args.input_module) else os.path.join(INPUT_DIR, args.input_module)
meta_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

print(f"📂 Using KO Data File:     {ko_file}")
print(f"📂 Using Module Data File: {mod_file}")
print(f"📂 Using Metadata File:    {meta_file}")
print(f"📂 Saving Outputs To:      {OUTPUT_DIR}")

# =========================================================
# 🔷 1. LOAD AND CLEAN METADATA, KO DATA & MODULE DATA
# =========================================================
print("Loading Metadata...")
if not os.path.exists(meta_file):
    print(f"❌ Error: Metadata file not found at {meta_file}")
    sys.exit(1)

meta = pd.read_csv(meta_file)
meta.columns = meta.columns.str.strip()

s_col = [c for c in meta.columns if "sample" in c.lower()][0] if any("sample" in c.lower() for c in meta.columns) else meta.columns[0]
h_col = [c for c in meta.columns if "host" in c.lower()][0] if any("host" in c.lower() for c in meta.columns) else meta.columns[1]

meta = meta.rename(columns={s_col: "SampleID", h_col: "Host"})
meta["SampleID"] = meta["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta["Host"] = meta["Host"].astype(str).str.strip()

# --- LOAD KO DATA FOR CHORD LINKS ---
print("Loading KO data for chord links...")
if not os.path.exists(ko_file):
    print(f"❌ Error: KO input file not found at {ko_file}")
    sys.exit(1)

try:
    df_ko = pd.read_csv(ko_file, sep="\t", header=None, names=["SampleID", "KO"])
    if df_ko.shape[1] < 2:
        df_ko = pd.read_csv(ko_file, sep=",", header=None, names=["SampleID", "KO"])
except Exception:
    df_ko = pd.read_csv(ko_file, sep=None, engine="python", header=None, names=["SampleID", "KO"])

df_ko["SampleID"] = df_ko["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
df_ko = df_ko.dropna(subset=["KO"])
df_ko = df_ko[df_ko["KO"] != "-"]
df_ko = df_ko.assign(KO=df_ko["KO"].astype(str).str.split(",")).explode("KO")
df_ko["KO"] = df_ko["KO"].str.replace(r"(?i)^ko:", "", regex=True).str.strip()

# --- LOAD MODULE DATA FOR INNER TRACKS ---
print("Loading Module data for inner tracks...")
if not os.path.exists(mod_file):
    print(f"❌ Error: Module input file not found at {mod_file}")
    sys.exit(1)

try:
    df_mod = pd.read_csv(mod_file, sep="\t", header=None, names=["SampleID", "Module"])
    if df_mod.shape[1] < 2:
        df_mod = pd.read_csv(mod_file, sep=",", header=None, names=["SampleID", "Module"])
except Exception:
    df_mod = pd.read_csv(mod_file, sep=None, engine="python", header=None, names=["SampleID", "Module"])

df_mod["SampleID"] = df_mod["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
df_mod = df_mod.dropna(subset=["Module"])
df_mod = df_mod[df_mod["Module"] != "-"]

# =========================================================
# 🔷 2. PROCESS KO OCCURRENCE & INTER-HOST SHARING MATRICES
# =========================================================
print("Processing KO shared matrices...")
ko_counts = df_ko.groupby(["KO", "SampleID"]).size().reset_index(name="Count")
df_host_ko = pd.merge(ko_counts, meta[["SampleID", "Host"]], on="SampleID", how="inner")

if df_host_ko.empty:
    print("❌ Error: Zero sample overlap between KO data and Metadata!")
    sys.exit(1)

df_summary_ko = df_host_ko.groupby(["Host", "KO"])["Count"].sum().reset_index()

ko_host_count = df_summary_ko.groupby("KO")["Host"].nunique().reset_index(name="n_host")
ko_host_dict = dict(zip(ko_host_count["KO"], ko_host_count["n_host"]))

mat = df_summary_ko.pivot(index="Host", columns="KO", values="Count").fillna(0)
mat_bin = (mat > 0).astype(int)

shared_matrix = mat_bin.dot(mat_bin.T)
np.fill_diagonal(shared_matrix.values, 0)

hosts = list(shared_matrix.index)
num_hosts = len(hosts)
ko_per_host = mat_bin.sum(axis=1).to_dict()

sectors = {host: max(ko_per_host.get(host, 1), 1) for host in hosts}

host_colors_list = sns.color_palette("Set2", n_colors=num_hosts).as_hex()
host_colors = dict(zip(hosts, host_colors_list))
sharing_colors = sns.color_palette("Spectral_r", n_colors=num_hosts).as_hex()

links = []
max_shared = shared_matrix.values.max() if shared_matrix.values.max() > 0 else 1

for i in range(num_hosts):
    for j in range(i + 1, num_hosts):
        h1, h2 = hosts[i], hosts[j]
        count = shared_matrix.loc[h1, h2]
        if count > 0:
            kos1 = set(mat_bin.columns[mat_bin.loc[h1] > 0])
            kos2 = set(mat_bin.columns[mat_bin.loc[h2] > 0])
            shared_kos = kos1.intersection(kos2)
            
            if shared_kos:
                avg_host = np.mean([ko_host_dict.get(k, 1) for k in shared_kos])
                color_index = int(np.round(avg_host)) - 1
                color_index = max(0, min(color_index, num_hosts - 1))
                link_color = sharing_colors[color_index]
            else:
                link_color = "#CCCCCC"

            links.append({
                "from": h1, "to": h2, "count": count,
                "color": link_color, "width": (count / max_shared) * 5
            })

# =========================================================
# 🔷 3. FETCH DETAILED KEGG MODULE HIERARCHY
# =========================================================
print("Fetching detailed KEGG module hierarchy...")
url = "https://rest.kegg.jp/get/br:ko00002"
try:
    response = requests.get(url, timeout=10)
    kegg_lines = response.text.splitlines()

    mod_ids, broad_cats, sub_cats = [], [], []
    current_broad, current_sub = "Unclassified", "Unclassified"

    for line in kegg_lines:
        if line.startswith("B"):
            current_broad = re.sub(r"^B\s*", "", line).strip()
            current_broad = re.sub(r"<[^>]+>", "", current_broad)
        elif line.startswith("C"):
            current_sub = re.sub(r"^C\s*", "", line).strip()
            current_sub = re.sub(r"<[^>]+>", "", current_sub)
        
        mod_match = re.search(r"\bM\d{5}\b", line)
        if mod_match:
            mod_ids.append(mod_match.group(0))
            broad_cats.append(current_broad)
            sub_cats.append(current_sub)

    df_cats = pd.DataFrame({
        "Module_ID": mod_ids,
        "Broad_Category": broad_cats,
        "Sub_Category": sub_cats
    }).drop_duplicates(subset=["Module_ID"], keep="first")
except Exception:
    print("⚠️ Web fetch failed. Falling back to default category tags...")
    df_cats = pd.DataFrame(columns=["Module_ID", "Broad_Category", "Sub_Category"])

mod_counts = df_mod.groupby(["SampleID", "Module"]).size().reset_index(name="Count")
mod_host = pd.merge(mod_counts, meta[["SampleID", "Host"]], on="SampleID", how="inner")
mod_summary = mod_host.groupby(["Host", "Module"])["Count"].sum().reset_index()

mod_summary = pd.merge(mod_summary, df_cats, left_on="Module", right_on="Module_ID", how="left")
mod_summary["Broad_Category"] = mod_summary["Broad_Category"].fillna("Unclassified").replace("", "Unclassified")
mod_summary["Sub_Category"] = mod_summary["Sub_Category"].fillna("Unclassified").replace("", "Unclassified")

host_max = mod_summary.groupby("Host")["Count"].transform("max")
host_max[host_max == 0] = 1.0
mod_summary["RelAbundance"] = mod_summary["Count"] / host_max

unique_broad = list(mod_summary["Broad_Category"].unique())
broad_palette = sns.color_palette("Set2", n_colors=max(1, len(unique_broad))).as_hex()
color_broad = dict(zip(unique_broad, broad_palette))
color_broad["Unclassified"] = "#E0E0E0"

unique_sub = list(mod_summary["Sub_Category"].unique())
sub_palette = sns.color_palette("Paired", n_colors=max(1, len(unique_sub))).as_hex()
color_sub = dict(zip(unique_sub, sub_palette))
color_sub["Unclassified"] = "#E0E0E0"

# =========================================================
# 🔷 4. RENDER CIRCOS PLOT
# =========================================================
print("Rendering Circos plot with clean separated layout...")
circos = Circos(sectors, space=3, start=0, end=360)

for sector in circos.sectors:
    host_name = sector.name
    
    # Outer Track: Plant Host Arcs (Radius 85-92)
    t_outer = sector.add_track((85, 92))
    t_outer.rect(sector.start, sector.end, fc=host_colors[host_name], ec="none")
    
    total_kos = ko_per_host.get(host_name, 0)
    interval = max(1, total_kos // 3)
    
    # Format tick marks clean and explicitly labeled
    t_outer.xticks_by_interval(
        interval=interval, 
        label_size=8, 
        label_orientation="vertical",
        label_formatter=lambda v: f"{int(v)} KOs"
    )
    
    # Sector Title Text placed clearly ABOVE the outer track (r=102)
    t_outer.text(
        f"{host_name}\n({total_kos} KOs)", 
        r=102, 
        size=10, 
        weight="bold",
        orientation="horizontal"
    )

    df_sub = mod_summary[mod_summary["Host"] == host_name].sort_values(
        by=["Broad_Category", "Sub_Category", "RelAbundance"],
        ascending=[True, True, False]
    )
    
    if len(df_sub) > 0:
        n_mod = len(df_sub)
        x_coords = np.linspace(sector.start, sector.end, n_mod + 1)
        
        # Middle Track: Broad Categories
        t_mid = sector.add_track((76, 83))
        for i, (_, row) in enumerate(df_sub.iterrows()):
            t_mid.rect(x_coords[i], x_coords[i+1], fc=color_broad[row["Broad_Category"]], ec="none")
            
        # Inner Track: Sub Categories
        t_inner = sector.add_track((67, 74))
        for i, (_, row) in enumerate(df_sub.iterrows()):
            t_inner.rect(x_coords[i], x_coords[i+1], fc=color_sub[row["Sub_Category"]], ec="none")

# Draw Chords
for link in links:
    sector_region1 = (link["from"], 0, sectors[link["from"]])
    sector_region2 = (link["to"], 0, sectors[link["to"]])
    circos.link(sector_region1, sector_region2, color=link["color"], alpha=0.4, r1=66, r2=66)

# =========================================================
# 🔷 5. BUILD PROPORTIONAL MULTI-COLUMN GRID
# =========================================================
num_sub = len(color_sub)
num_broad = len(color_broad)

# Dynamic canvas height scaling to handle long subcategory lists
dynamic_height = max(16, 8 + (num_sub // 4) * 0.45)
fig = plt.figure(figsize=(26, dynamic_height))

# Split Top (Circos Plot) and Bottom (Legends Panel)
gs_main = gridspec.GridSpec(2, 1, height_ratios=[1.8, 1.0], hspace=0.20)

# 1. Top Subplot: Circos Polar Plot
ax_circos = fig.add_subplot(gs_main[0], projection="polar")
circos.plotfig(ax=ax_circos)

# 2. Bottom Grid: Proportional width ratios shift Broad Categories left
# [0.20 (Hosts/Sharing), 0.25 (Broad Cats), 0.55 (Sub Cats)]
gs_bottom = gridspec.GridSpecFromSubplotSpec(
    1, 3, 
    subplot_spec=gs_main[1], 
    width_ratios=[0.20, 0.25, 0.55], 
    wspace=0.15
)

ax_leg1 = fig.add_subplot(gs_bottom[0, 0])  # Hosts & Sharing
ax_leg2 = fig.add_subplot(gs_bottom[0, 1])  # Broad Categories
ax_leg3 = fig.add_subplot(gs_bottom[0, 2])  # Sub Categories

ax_leg1.axis("off")
ax_leg2.axis("off")
ax_leg3.axis("off")

# ---------------------------------------------------------
# LEGEND 1: Plant Hosts & KO Sharing (Compact Left)
# ---------------------------------------------------------
legend_hosts = [Line2D([0], [0], color='none', label=r'$\bf{Plant\ Hosts}$')]
for h, c in host_colors.items():
    legend_hosts.append(mpatches.Patch(color=c, label=h))

leg_host_obj = ax_leg1.legend(handles=legend_hosts, loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=9.5, frameon=False)
ax_leg1.add_artist(leg_host_obj)

legend_sharing = [Line2D([0], [0], color='none', label=r'$\bf{KO\ Sharing}$')]
for i in range(num_hosts):
    lbl = f"{i+1} Host" if i == 0 else f"{i+1} Hosts"
    legend_sharing.append(mpatches.Patch(color=sharing_colors[i], label=lbl))

sharing_anchor_y = max(0.1, 1.0 - (len(host_colors) * 0.08))
leg_sharing_obj = ax_leg1.legend(handles=legend_sharing, loc="upper left", bbox_to_anchor=(0, sharing_anchor_y), fontsize=9.0, frameon=False)
ax_leg1.add_artist(leg_sharing_obj)

# ---------------------------------------------------------
# LEGEND 2: Broad Categories (Shifted Left)
# ---------------------------------------------------------
legend_broad = [Line2D([0], [0], color='none', label=r'$\bf{Broad\ Categories}$')]
for bcat, c in color_broad.items():
    label_short = bcat[:32] + "..." if len(bcat) > 34 else bcat
    legend_broad.append(mpatches.Patch(color=c, label=label_short))

ax_leg2.legend(handles=legend_broad, loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=8.5, frameon=False, ncol=1)

# ---------------------------------------------------------
# LEGEND 3: Sub Categories (4 Columns for Horizontal Fit)
# ---------------------------------------------------------
legend_sub = [Line2D([0], [0], color='none', label=r'$\bf{Sub\ Categories}$')]
for scat, c in color_sub.items():
    label_short = scat[:26] + "..." if len(scat) > 28 else scat
    legend_sub.append(mpatches.Patch(color=c, label=label_short))

ncol_sub = 4 if num_sub >= 30 else 3
ax_leg3.legend(handles=legend_sub, loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=8.0, frameon=False, ncol=ncol_sub)

png_output_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_KO_circos_plot.png")
plt.savefig(png_output_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ PNG Saved: {png_output_path}")

# =========================================================
# 🔷 6. SAVE INTERACTIVE HTML MATRIX SUMMARY
# =========================================================
print("Generating interactive HTML summary...")
fig_html = px.imshow(
    shared_matrix,
    labels=dict(x="Host", y="Host", color="Shared KO Count"),
    x=shared_matrix.columns, y=shared_matrix.index,
    color_continuous_scale="Spectral_r",
    title="Shared KO Numbers Between Plant Hosts (Interactive Summary)"
)
fig_html.update_layout(
    xaxis=dict(tickangle=-45, tickfont=dict(size=11, color="black")),
    yaxis=dict(tickfont=dict(size=11, color="black")),
    template="plotly_white"
)

html_output_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_KO_shared_matrix.html")
fig_html.write_html(html_output_path)
print(f"✅ HTML Saved: {html_output_path}")

print("🎉 Completed successfully!")
