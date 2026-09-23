#!/usr/bin/env python3
"""
EC Circos Diagram and Feature Sharing Pipeline (Proportional Legend Grid)
-------------------------------------------------------------------------
Resolves layout gaps by using custom width ratios for the bottom legend
panel, shifting Broad Categories leftward and balancing space across columns.
"""

import os
import re
import sys
import argparse
import requests
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

import seaborn as sns
from pycirclize import Circos
import plotly.express as px

# =========================================================
# 1. COMMAND LINE ARGUMENTS
# =========================================================
parser = argparse.ArgumentParser(
    description="EC Circos Diagram Pipeline (Proportional Layout)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("-i", "--input", type=str, default="all_sample_ec_MJ.tsv", help="EC long-format file")
parser.add_argument("-m", "--metadata", type=str, default="Plant-host2.csv", help="Metadata CSV")
parser.add_argument("-o", "--output-dir", type=str, default="/Metrico_project/output", help="Output directory")
parser.add_argument("--prefix", type=str, default="ec_circos", help="Prefix for output files")

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY & FILE PATHS
# =========================================================
BASE_DIR = "/Metrico_project"
INPUT_DIR = os.path.join(BASE_DIR, "input_file")
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

ec_file = args.input if os.path.isabs(args.input) else os.path.join(INPUT_DIR, args.input)
meta_file = args.metadata if os.path.isabs(args.metadata) else os.path.join(INPUT_DIR, args.metadata)

# =========================================================
# 🔷 1. LOAD & CLEAN DATA
# =========================================================
if not os.path.exists(ec_file) or not os.path.exists(meta_file):
    print("❌ Error: Input files not found!")
    sys.exit(1)

try:
    df = pd.read_csv(ec_file, sep="\t", header=None, names=["SampleID", "EC"])
    if df.shape[1] < 2:
        df = pd.read_csv(ec_file, sep=",", header=None, names=["SampleID", "EC"])
except Exception:
    df = pd.read_csv(ec_file, sep=None, engine="python", header=None, names=["SampleID", "EC"])

meta = pd.read_csv(meta_file)
meta.columns = meta.columns.str.strip()

s_col = [c for c in meta.columns if "sample" in c.lower()][0] if any("sample" in c.lower() for c in meta.columns) else meta.columns[0]
h_col = [c for c in meta.columns if "host" in c.lower()][0] if any("host" in c.lower() for c in meta.columns) else meta.columns[1]

meta = meta.rename(columns={s_col: "SampleID", h_col: "Host"})
meta["SampleID"] = meta["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
meta["Host"] = meta["Host"].astype(str).str.strip()

df["SampleID"] = df["SampleID"].astype(str).str.strip().str.replace(r"_eggnog.*", "", regex=True)
df = df.dropna(subset=["EC"])
df = df[df["EC"] != "-"]
df = df.assign(EC=df["EC"].astype(str).str.split(",")).explode("EC")
df["EC"] = df["EC"].str.replace(r"(?i)^ec:", "", regex=True).str.strip()

# =========================================================
# 🔷 2. BUILD OCCURRENCE & SHARING MATRICES
# =========================================================
ec_counts = df.groupby(["EC", "SampleID"]).size().reset_index(name="Count")
df_host = pd.merge(ec_counts, meta[["SampleID", "Host"]], on="SampleID", how="inner")

df_summary = df_host.groupby(["Host", "EC"])["Count"].sum().reset_index()

ec_host_count = df_summary.groupby("EC")["Host"].nunique().reset_index(name="n_host")
ec_host_dict = dict(zip(ec_host_count["EC"], ec_host_count["n_host"]))

mat = df_summary.pivot(index="Host", columns="EC", values="Count").fillna(0)
mat_bin = (mat > 0).astype(int)

shared_matrix = mat_bin.dot(mat_bin.T)
np.fill_diagonal(shared_matrix.values, 0)

hosts = list(shared_matrix.index)
num_hosts = len(hosts)
ec_per_host = mat_bin.sum(axis=1).to_dict()

sectors = {host: max(ec_per_host.get(host, 1), 1) for host in hosts}

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
            ecs1 = set(mat_bin.columns[mat_bin.loc[h1] > 0])
            ecs2 = set(mat_bin.columns[mat_bin.loc[h2] > 0])
            shared_ecs = ecs1.intersection(ecs2)
            
            if shared_ecs:
                avg_host = np.mean([ec_host_dict.get(ec, 1) for ec in shared_ecs])
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
# 🔷 3. FETCH KEGG HIERARCHY
# =========================================================
url = "https://rest.kegg.jp/get/br:ko01000"
try:
    response = requests.get(url, timeout=10)
    kegg_lines = response.text.splitlines()

    ec_ids, broad_cats, sub_cats = [], [], []
    current_broad, current_sub = "", ""

    for line in kegg_lines:
        if line.startswith("B"):
            raw_broad = re.sub(r"<[^>]+>", "", line[1:]).strip()
            current_broad = re.sub(r"^[\d\.]+\s*", "", raw_broad)
        elif line.startswith("C"):
            raw_sub = re.sub(r"<[^>]+>", "", line[1:]).strip()
            current_sub = re.sub(r"^[\d\.]+\s*", "", raw_sub)
        
        ec_match = re.search(r"\b\d+\.\d+\.\d+\.[\d\-]+\b", line)
        if ec_match:
            ec_ids.append(ec_match.group(0))
            broad_cats.append(current_broad)
            sub_cats.append(current_sub)

    df_cats = pd.DataFrame({
        "EC_ID": ec_ids, "Broad_Category": broad_cats, "Sub_Category": sub_cats
    }).drop_duplicates(subset=["EC_ID"], keep="first")
except Exception:
    df_cats = pd.DataFrame(columns=["EC_ID", "Broad_Category", "Sub_Category"])

ec_summary = pd.merge(df_summary, df_cats, left_on="EC", right_on="EC_ID", how="left")
ec_summary["Broad_Category"] = ec_summary["Broad_Category"].fillna("Unclassified").replace("", "Unclassified")
ec_summary["Sub_Category"] = ec_summary["Sub_Category"].fillna("Unclassified").replace("", "Unclassified")

host_max = ec_summary.groupby("Host")["Count"].transform("max")
host_max[host_max == 0] = 1.0
ec_summary["RelAbundance"] = ec_summary["Count"] / host_max

unique_broad = list(ec_summary["Broad_Category"].unique())
broad_palette = sns.color_palette("Set2", n_colors=max(1, len(unique_broad))).as_hex()
color_broad = dict(zip(unique_broad, broad_palette))
color_broad["Unclassified"] = "#E0E0E0"

unique_sub = list(ec_summary["Sub_Category"].unique())
sub_palette = sns.color_palette("Paired", n_colors=max(1, len(unique_sub))).as_hex()
color_sub = dict(zip(unique_sub, sub_palette))
color_sub["Unclassified"] = "#E0E0E0"

# =========================================================
# 🔷 4. RENDER CIRCOS PLOT
# =========================================================
print("Rendering Circos plot with adjusted proportional layout...")
circos = Circos(sectors, space=3, start=0, end=360)

for sector in circos.sectors:
    host_name = sector.name
    
    # Outer Track: Plant Host Arcs
    t_outer = sector.add_track((85, 92))
    t_outer.rect(sector.start, sector.end, fc=host_colors[host_name], ec="none")
    
    total_ecs = ec_per_host.get(host_name, 0)
    interval = max(1, total_ecs // 3)
    
    t_outer.xticks_by_interval(
        interval=interval, 
        label_size=8, 
        label_orientation="vertical",
        label_formatter=lambda v: f"{int(v)} ECs"
    )
    
    t_outer.text(
        f"{host_name}\n({total_ecs} ECs)", 
        r=102, 
        size=10, 
        weight="bold",
        orientation="horizontal"
    )

    df_sub = ec_summary[ec_summary["Host"] == host_name].sort_values(
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

dynamic_height = max(16, 8 + (num_sub // 4) * 0.45)
fig = plt.figure(figsize=(26, dynamic_height))

# Split Top (Circos Plot) and Bottom (Legends Panel)
gs_main = gridspec.GridSpec(2, 1, height_ratios=[1.8, 1.0], hspace=0.20)

# 1. Top Subplot: Circos Polar Plot
ax_circos = fig.add_subplot(gs_main[0], projection="polar")
circos.plotfig(ax=ax_circos)

# 2. Bottom Grid: Custom width ratios shift Broad Categories left
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
# LEGEND 1: Plant Hosts & EC Sharing (Compact Left)
# ---------------------------------------------------------
legend_hosts = [Line2D([0], [0], color='none', label=r'$\bf{Plant\ Hosts}$')]
for h, c in host_colors.items():
    legend_hosts.append(mpatches.Patch(color=c, label=h))

leg_host_obj = ax_leg1.legend(handles=legend_hosts, loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=9.5, frameon=False)
ax_leg1.add_artist(leg_host_obj)

legend_sharing = [Line2D([0], [0], color='none', label=r'$\bf{EC\ Sharing}$')]
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
# LEGEND 3: Sub Categories (4 Columns to maximize horizontal fit)
# ---------------------------------------------------------
legend_sub = [Line2D([0], [0], color='none', label=r'$\bf{Sub\ Categories}$')]
for scat, c in color_sub.items():
    label_short = scat[:26] + "..." if len(scat) > 28 else scat
    legend_sub.append(mpatches.Patch(color=c, label=label_short))

ncol_sub = 4 if num_sub >= 30 else 3
ax_leg3.legend(handles=legend_sub, loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=8.0, frameon=False, ncol=ncol_sub)

plt.savefig(os.path.join(OUTPUT_DIR, f"{args.prefix}_EC_circos_plot.png"), dpi=300, bbox_inches="tight")
plt.close()

# =========================================================
# 🔷 6. SAVE INTERACTIVE HTML MATRIX SUMMARY
# =========================================================
fig_html = px.imshow(
    shared_matrix,
    labels=dict(x="Host", y="Host", color="Shared EC Count"),
    x=shared_matrix.columns, y=shared_matrix.index,
    color_continuous_scale="Spectral_r",
    title="Shared EC Numbers Between Plant Hosts (Interactive Summary)"
)
fig_html.update_layout(
    xaxis=dict(tickangle=-45, tickfont=dict(size=11, color="black")),
    yaxis=dict(tickfont=dict(size=11, color="black")),
    template="plotly_white"
)

html_output_path = os.path.join(OUTPUT_DIR, f"{args.prefix}_shared_matrix.html")
fig_html.write_html(html_output_path)

print("🎉 Processed proportional, balanced layout successfully!")
