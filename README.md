# 🧬 Metrico_project: Metagenomic Functional Visualization Pipeline

Welcome to **Metrico_project**! This toolkit provides a streamlined Python and Docker pipeline for processing and visualizing metagenomic functional annotations across plant hosts and environmental samples (including CAZy Classes, GH Families, COG Categories, EC Matrices, KO Matrices, EC Circos Plots, KO Circos Plots, MetaCyc Pathway Bubble Plots, and PCoA/PERMANOVA Ordination).

---

## ⚡ Quick Start: Container & Repository Setup

Before running any visualization modules, pull the pre-packaged Docker image from Docker Hub and clone the repository:

### 1. Pull the Docker Image from Docker Hub
No manual building is required. Pull the ready-to-run container directly:
```bash
docker pull mantrilabnabi/metrico_viz:v1
```

### 2. Clone the Metrico Repository
Large input matrices in input_file/ are tracked using Git LFS. You must install and initialize git-lfs before cloning or run git lfs pull after cloning to download actual data files instead of pointer files.
```bash
# Install and initialize Git LFS (Ubuntu/Debian)
sudo apt-get update && sudo apt-get install -y git-lfs
git lfs install

# Clone repository and navigate into folder
git clone [https://github.com/cbl-nabiadmin/Metrico_project.git](https://github.com/cbl-nabiadmin/Metrico_project.git)
cd Metrico_project

# Pull actual Git LFS binary datasets
git lfs pull
```

### 3. Quick Test Run (Using Pre-packaged Data)
Execute the master runner immediately with default test datasets to verify your setup:
```bash
sudo docker run --rm \
  -v $(pwd):/Metrico_project \
  mantrilabnabi/metrico_viz:v1 python3 /Metrico_project/run_all.py
```

> **Note on Prerequisites:** If you plan to analyze your own custom samples rather than the provided example data, you must first convert your raw annotation outputs (from eggNOG-mapper or HUMAnN) into compatible matrix files using the pre-processing scripts detailed in [⚙️ Prerequisites & Data Pre-Processing Tutorial](#️-prerequisites--data-pre-processing-tutorial).

---

## 📌 Navigation Index

* [⚡ Quick Start: Container & Repository Setup](#-quick-start-container--repository-setup)
* [🧭 Workflow Fast-Track & Intermediate Entry Points](#-workflow-fast-track--intermediate-entry-points)
* [⚙️ Prerequisites & Data Pre-Processing Tutorial](#️-prerequisites--data-pre-processing-tutorial)
* [📁 Project Directory Structure](#-project-directory-structure)
* [🛠️ Step 1: Unpack & Load the Docker Container](#️-step-1-unpack--load-the-docker-container)
* [🚀 Step 2: Running the Master Pipeline (`run_all.py`)](#-step-2-running-the-master-pipeline-run_allpy)
* [📋 Full Parameter Reference Table for `run_all.py`](#-full-parameter-reference-table-for-run_allpy)
* [🎯 Running Individual Scripts Separately](#-running-individual-scripts-separately)
* [📊 Output File Deliverables](#-output-file-deliverables)

---

## 🧭 Workflow Fast-Track & Intermediate Entry Points

You do not need to run every pre-processing step or tool if you already have processed matrices or only want specific plots. Use the table below to determine where to start:

| Your Current Starting Point | Can You Skip Pre-Processing? | Actions Required |
| :--- | :---: | :--- |
| **I want to test the pipeline using provided example data** | **YES** *(Skip Tasks 1–5)* | Pull the Docker image, clone the repo, and run `run_all.py` directly. |
| **I already have formatted matrices** *(e.g., `all_ko_matrix.tsv`, `cazy_class_matrix.tsv`)* | **YES** *(Skip Tasks 1–5)* | Place your files into the `input_file/` directory and jump directly to [Step 2](#-step-2-running-the-master-pipeline-run_allpy) or [Individual Scripts](#-running-individual-scripts-separately). |
| **I only want to run specific plots** *(e.g., only CAZy or COG)* | **PARTIAL** | Only execute the specific pre-processing tasks required for your target plot, then run the corresponding Individual Script. |
| **I have raw eggNOG (`.emapper.annotations`) or HUMAnN (`pathabundance.tsv`) outputs** | **NO** | Complete the relevant pre-processing tasks below before running any visualization scripts. |

---

## ⚙️ Prerequisites & Data Pre-Processing Tutorial

> **⚠️ CRITICAL REQUIREMENT:** Unless you already have pre-formatted matrices, you must first process your raw FASTQ files using HUMAnN and your assembled FASTA files using eggNOG-mapper (following assembly with SPAdes). Once processed, convert the resulting outputs into the required matrix formats using the commands below.

### Pre-Processing Workflow Map

| Task | Input File | Target Visualization | Required Output Name |
| :--- | :--- | :--- | :--- |
| **Combine Annotations** | `*.emapper.annotations` | Central Master Dataset | `all_with_source.annotations` |
| **Circos KO & Modules** | `all_with_source.annotations` | Circular KO/Module Plot | `all_sample_ko.tsv`, `all_module.tsv` |
| **Circos EC Numbers** | `all_with_source.annotations` | Circular EC Chord Diagram | `all_sample_ec.tsv` |
| **HUMAnN Abundance** | `*_pathabundance.tsv` | Pathway Driver Bubble Plot | `all-pathabundance_cpm.tsv` |
| **COG Category Matrix** | `all_with_source.annotations` | COG Functional Barplot | `all_cog_matrix.csv` |
| **KO Matrix** | `all_with_source.annotations` | Functional KO Heatmap | `all_ko_matrix.tsv` |
| **EC Matrix** | `all_with_source.annotations` | Enzyme EC Heatmap | `all_ec_matrix.tsv` |
| **CAZy Profile Matrix** | `all_with_source.annotations` | CAZy Class/Family Heatmap | `cazy_class_matrix.tsv`, `gh_family_filtered.csv` / `gt_family_filtered.csv` |

---

### Step-by-Step Pre-Processing Commands

#### Task 1: Combine eggNOG Annotations
Merge individual sample subfolder `.emapper.annotations` files into a single master annotation file with sample labels attached. Save this script as `combined_results.sh` and run it within the directory containing sample subdirectories.

```bash
#!/bin/bash
OUT="all_with_source.annotations"
DIR="/Path/to/emapper-result/samplewise"

# Get first annotation file to extract header
first_file=$(find "$DIR" -type f -name "*.emapper.annotations" | head -n 1)

# Extract header and add Sample column
echo -e "Sample\t$(grep "^#query" "$first_file" | sed 's/^#//')" > "$OUT"

# Loop through each subdirectory
for d in "$DIR"/*/; do
    sample=$(basename "$d")
    f=$(find "$d" -maxdepth 1 -name "*.emapper.annotations")
    [ -z "$f" ] && continue
    grep -v "^#" "$f" | awk -v s="$sample" '{print s"\t"$0}'
done >> "$OUT"
```

**Run:**
```bash
bash combined_results.sh
```
**Output:** `all_with_source.annotations`

---

#### Task 2: Circos Inputs (KO & Modules, EC Numbers)
Extract KO IDs, unroll comma-separated KEGG Modules into long-format files, and extract annotated EC numbers. Save this script as `circos_script.sh` and run it in the directory containing `all_with_source.annotations`.

```bash
#!/bin/bash
# 1. Extract KO IDs
awk -F '\t' '{print $1"\t"$13}' all_with_source.annotations > all_sample_ko.tsv

# 2. Extract & Unroll KEGG Modules
awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ && $15 != "-" {
    split($15, mod, ",");
    for(i in mod){
        gsub(/^ +| +$/, "", mod[i]);
        split($1, a, "_");
        print a[1], mod[i]
    }
}' all_with_source.annotations > all_module.tsv

# 3. Extract sample IDs and EC numbers
awk -F '\t' '$12 != "" && $12 != "-" {print $1"\t"$12}' all_with_source.annotations > all_sample_ec.tsv
```

**Run:**
```bash
bash circos_script.sh
```
**Output:** `all_sample_ko.tsv`, `all_module.tsv`, `all_sample_ec.tsv`

---

#### Task 3: HUMAnN Pathway CPM Normalization
Merge per-sample HUMAnN pathway abundances and normalize to Counts Per Million (CPM). Save this script as `humann_rejoin.sh`.

```bash
#!/bin/bash
BASE="/path/to/human-result"
JOIN_DIR="/tmp/staged_pathabundance"
OUTDIR="./joined_output"
mkdir -p "$JOIN_DIR" "$OUTDIR"

# Stage files into staging folder
for FILE in "$BASE"/*/*_merged_pathabundance.tsv; do
    [ -f "$FILE" ] && cp "$FILE" "$JOIN_DIR/"
done

# Merge and normalize tables
humann_join_tables -i "$JOIN_DIR" -o "$OUTDIR/all-pathabundance.tsv" --file_name pathabundance
humann_renorm_table -i "$OUTDIR/all-pathabundance.tsv" -o "$OUTDIR/all-pathabundance_cpm.tsv" --units cpm --update-snames

# Cleanup staging directory
rm -rf "$JOIN_DIR"
```

**Run:**
```bash
bash humann_rejoin.sh
```
**Output:** `all-pathabundance_cpm.tsv`

---

#### Task 4: Build COG, EC & KO Matrices
Parse primary COG categories, KO IDs, EC numbers, and CAZy annotations into wide-format matrices. Save as `matrix_build.sh`.

```bash
#!/bin/bash
# 1. Extract COG categories and pivot
awk -F '\t' '$8 != "-" && $8 != "" { print $1 "\t" substr($8,1,1) }' all_with_source.annotations > all_cog_origin2.tsv
sort all_cog_origin2.tsv | uniq -c | awk '{print $2"\t"$3"\t"$1}' > all_cog_counts2.tsv

awk '{
    sample=$1; cog=$2; val=$3;
    count[sample"\t"cog]=val; samples[sample]=1; cogs[cog]=1;
}
END {
    printf "SampleID";
    for (c in cogs) printf "\t" c;
    printf "\n";
    for (s in samples) {
        printf s;
        for (c in cogs) { key=s"\t"c; printf "\t%d", (key in count ? count[key] : 0); }
        printf "\n";
    }
}' all_cog_counts2.tsv > all_cog_matrix.csv

# 2. Extract KO IDs and pivot
awk -F '\t' '$13 != "-" {
    split($13, a, ",");
    for (i in a) { gsub("ko:", "", a[i]); print $1 "\t" a[i]; }
}' all_with_source.annotations > all-ko-origin.tsv

sort all-ko-origin.tsv | uniq -c | awk '{print $2"\t"$3"\t"$1}' > all_ko_counts.tsv

awk '{ key=$2"\t"$1; count[key]=$3; kos[$2]=1; samples[$1]=1; }
END {
    printf "KO"; for (s in samples) printf "\t"s; printf "\n";
    for (k in kos) {
        printf k;
        for (s in samples) { key=k"\t"s; printf "\t%d", (key in count ? count[key] : 0); }
        printf "\n";
    }
}' all_ko_counts.tsv > all_ko_matrix.tsv

# 3. Extract EC numbers and pivot
awk -F '\t' '$12 != "-" && $12 != "" {
    split($12, a, ",");
    for (i in a) { gsub(/^ +| +$/, "", a[i]); if (a[i] != "") print $1 "\t" a[i]; }
}' all_with_source.annotations > all_ec_origin.tsv

sort all_ec_origin.tsv | uniq -c | awk '{print $2"\t"$3"\t"$1}' > all_ec_counts.tsv

awk '{ key=$2"\t"$1; count[key]=$3; ecs[$2]=1; samples[$1]=1; }
END {
    printf "EC"; for (s in samples) printf "\t"s; printf "\n";
    for (e in ecs) {
        printf e;
        for (s in samples) { key=e"\t"s; printf "\t%d", (key in count ? count[key] : 0); }
        printf "\n";
    }
}' all_ec_counts.tsv > all_ec_matrix.tsv

# 4. Extract CAZy annotations
awk -F '\t' 'NR>1 { if($20 != "-" && $20 != "") print $1 "\t" $20 }' all_with_source.annotations | \
awk -F '\t' '{
    n=split($2, a, ",");
    for(i=1;i<=n;i++){ gsub(/^ +| +$/, "", a[i]); if(a[i] != "" && a[i] != "-") print $1 "\t" a[i]; }
}' > all_cazy_expanded.tsv

awk -F '\t' '{ class=$2; sub(/[0-9]+.*/, "", class); print $1 "\t" $2 "\t" class; }' all_cazy_expanded.tsv > all_cazy_class.tsv
export LC_ALL=C
sort all_cazy_class.tsv | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > all_cazy_count.csv
```

**Run:**
```bash
bash matrix_build.sh
```
**Output:** `all_cog_matrix.csv`, `all_ko_matrix.tsv`, `all_ec_matrix.tsv`, `all_cazy_count.csv`

---

#### Task 5: CAZy Profile & Subfamily Matrices
Convert long-format CAZy counts to a wide matrix using Python (`cazy_matrix_convert.py`):

```python
#!/usr/bin/env python3
import pandas as pd
import sys

input_file = "all_cazy_count.csv"
output_file = "cazy_class_matrix.tsv"

try:
    df = pd.read_csv(input_file, sep=None, engine="python", header=None)
    if str(df.iloc[0, 0]).lower().startswith("sample"):
        df = df.iloc[1:]
    df.columns = ["SampleID", "Family", "Class", "Count"]
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)

df["SampleID"] = df["SampleID"].astype(str).str.strip()
df["Class"] = df["Class"].astype(str).str.strip()
df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0)

matrix_df = df.groupby(["SampleID", "Class"])["Count"].sum().unstack(fill_value=0)
matrix_df.index.name = "Sample"
matrix_df.to_csv(output_file, sep="\t")

print(f"✅ Conversion complete! Matrix saved to: {output_file}")
```

**Run:**
```bash
python3 cazy_matrix_convert.py
```
**Output:** `cazy_class_matrix.tsv`

To extract specific CAZy subfamilies (e.g., GH, GT, CBM), save and run `filter_cazy_class.py`:

```python
#!/usr/bin/env python3
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Generic CAZy Class Extractor")
parser.add_argument("-i", "--input", type=str, required=True, help="Input raw CAZy CSV/TSV file")
parser.add_argument("-c", "--cazy-class", type=str, default="GH", help="Target CAZy class (e.g., GH, GT)")
parser.add_argument("-o", "--output", type=str, default=None, help="Output CSV path")

args = parser.parse_args()
target_class = args.cazy_class.upper().strip()
output_file = args.output if args.output else f"{target_class.lower()}_family_filtered.csv"

df = pd.read_csv(args.input, sep=None, engine="python", header=None)
if str(df.iloc[0, 0]).lower().startswith("sample"):
    df = df.iloc[1:]
df = df.iloc[:, :4]
df.columns = ["SampleID", "Family", "Class", "Count"]

df["SampleID"] = df["SampleID"].astype(str).str.strip()
df["Family"] = df["Family"].astype(str).str.strip()
df["Class"] = df["Class"].astype(str).str.upper().str.strip()
df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0).astype(int)

filtered_df = df[(df["Class"] == target_class) | (df["Family"].str.startswith(target_class))].copy()
family_col_name = f"{target_class}_Family"
filtered_df = filtered_df.rename(columns={"Family": family_col_name})
filtered_df = filtered_df[["Count", "SampleID", family_col_name]]
filtered_df.to_csv(output_file, index=False)

print(f"✅ Extracted {len(filtered_df)} entries for '{target_class}'. Saved to: {output_file}")
```

**Run:**
```bash
python3 filter_cazy_class.py -i all_cazy_count.csv -c GH -o gh_family_filtered.csv
```
**Output:** `gh_family_filtered.csv`

---

## 📁 Project Directory Structure

Once pre-processing is complete, place your generated files into the `input_file/` directory as structured below:

```text
Metrico_project/
├── input_file/                        # Store input matrices & metadata here
│   ├── cazy_class_matrix.tsv          # Example CAZy Class matrix
│   ├── gh_family_filtered.csv         # Example GH Family count matrix
│   ├── all_cog_matrix.csv             # Example COG functional matrix
│   ├── all_ec_matrix.tsv              # Example EC count matrix
│   ├── all_sample_ec.tsv              # Example long-format EC file
│   ├── all_ko_matrix.tsv              # Example KO count matrix
│   ├── all_sample_ko.tsv              # Example long-format KO file
│   ├── all-pathabundance_cpm.tsv      # Example MetaCyc pathway abundance file
│   ├── MetaCyc-pathways.txt           # MetaCyc ontology mapping file
│   └── Plant-host.csv                 # Example plant host metadata
├── output/                            # Output results directory
├── cazy_heatmap.py                    # CAZy Class analysis script
├── gh_family_heatmap.py               # GH Family analysis script
├── cog_abundance_plot.py              # COG analysis script
├── ec_heatmap.py                      # EC Matrix analysis script
├── ec_circos_plot.py                  # EC Circos analysis script
├── ko_heatmap.py                      # KO Matrix analysis script
├── ko_circos_plot.py                  # KO Circos analysis script
├── pathway_bubble_plot.py             # Pathway Bubble plot script
├── ordination_permanova.py           # PCoA & PERMANOVA ordination script
└── run_all.py                         # Master pipeline wrapper script
```

---

## 🛠️ Step 1: Unpack & Load the Docker Container

The pipeline environment can be pulled directly from Docker Hub or loaded from a compressed archive (`.tar.gz`).

### Option A: Pull Pre-built Container from Docker Hub (Recommended)
```bash
docker pull mantrilabnabi/metrico_viz:v1
```

### Option B: Load from Local Archive File
```bash
cd /path/to/image/archive
sudo docker load -i metrico_viz:v1.tar.gz
```

Verify that the image has loaded successfully:
```bash
sudo docker images
```
*(You should see `mantrilabnabi/metrico_viz:v1` or `metrico_viz:v1` listed under REPOSITORY).*

---

## 🚀 Step 2: Running the Master Pipeline (`run_all.py`)

### 1. Default Run (Using Example Data)
Example input datasets are pre-packaged inside the `input_file/` directory. Running `run_all.py` without extra flags processes all default example datasets sequentially and outputs figures into `output/`:

```bash
sudo docker run --rm \
  -v $(pwd):/Metrico_project \
  mantrilabnabi/metrico_viz:v1 python3 /Metrico_project/run_all.py
```

### 2. Full Master Command (Custom Data Across All 9 Modules)
Pass custom input files directly via command-line arguments:

```bash
sudo docker run --rm \
  -v $(pwd):/Metrico_project \
  mantrilabnabi/metrico_viz:v1 python3 /Metrico_project/run_all.py \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output \
    --prefix Metrico \
    --mod /Metrico_project/input_file/all_module.tsv \
    --cazy-class-input /Metrico_project/input_file/cazy_class_matrix.tsv \
    --gh-family-input /Metrico_project/input_file/gh_family_filtered.csv \
    --cog-input /Metrico_project/input_file/all_cog_matrix.csv \
    --ec-matrix-input /Metrico_project/input_file/all_ec_matrix.tsv \
    --ec-circos-input /Metrico_project/input_file/all_sample_ec.tsv \
    --ko-matrix-input /Metrico_project/input_file/all_ko_matrix.tsv \
    --ko-circos-input /Metrico_project/input_file/all_sample_ko.tsv \
    --pathway-input /Metrico_project/input_file/all-pathabundance_cpm.tsv \
    --ordination-input /Metrico_project/input_file/all_ko_matrix.tsv \
    --ordination-metric aitchison \
    --ellipse-min-n 3 \
    --permutations 999
```

---

## 📋 Full Parameter Reference Table for `run_all.py`

| Parameter Flag | Description | Default Value |
| :--- | :--- | :--- |
| `--rm` | Container Cleanup Flag: Removes anonymous container filesystem upon exit | Docker Runtime Flag |
| `-m` / `--metadata` | Common metadata CSV file mapping samples to experimental groups | `Plant-host.csv` |
| `-o` / `--output-dir` | Target directory location for all generated output files | `/Metrico_project/output` |
| `--prefix` | Global prefix tag prepended to saved output filenames | `Metrico` |
| `--mod` / `--input-module` | KEGG Module file required for KO Circos plot | `all_module.tsv` |
| `--cazy-class-input` | Input file for CAZy Class heatmap (`cazy_heatmap.py`) | `cazy_class_matrix.tsv` |
| `--gh-family-input` | Input file for GH Family heatmap (`gh_family_heatmap.py`) | `gh_family_filtered.csv` |
| `--cog-input` | Input file for COG distribution (`cog_abundance_plot.py`) | `all_cog_matrix.csv` |
| `--ec-matrix-input` | Input file for EC heatmap (`ec_heatmap.py`) | `all_ec_matrix.tsv` |
| `--ec-circos-input` | Input file for EC Circos plot (`ec_circos_plot.py`) | `all_sample_ec.tsv` |
| `--ko-matrix-input` | Input file for KO heatmap (`ko_heatmap.py`) | `all_ko_matrix.tsv` |
| `--ko-circos-input` | Input file for KO Circos plot (`ko_circos_plot.py`) | `all_sample_ko.tsv` |
| `--pathway-input` | Input file for MetaCyc bubble plot (`pathway_bubble_plot.py`) | `all-pathabundance_cpm.tsv` |
| `--ordination-input` | Input file for PCoA/PERMANOVA (`ordination_permanova.py`) | `all_ko_matrix.tsv` |
| `--ordination-metric` | Distance metric for ordination (`aitchison` or `bray`) | `aitchison` |
| `--ellipse-min-n` | Minimum group size required to draw 95% confidence ellipses | `3` |
| `--permutations` | Number of Monte Carlo permutations for PERMANOVA testing | `999` |

---

## 🎯 Running Individual Scripts Separately

You can execute any analysis script independently using `-i` (input), `-m` (metadata), `-o` (output directory), and `--prefix` (filename prefix):

### A. CAZy Class Heatmap
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/cazy_heatmap.py \
    -i /Metrico_project/input_file/cazy_class_matrix.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output --prefix CAZy
```

### B. GH Family Heatmap
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/gh_family_heatmap.py \
    -i /Metrico_project/input_file/gh_family_filtered.csv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -n 30 -o /Metrico_project/output --prefix GH_Family
```

### C. COG Categories Bar Plots & Heatmap
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/cog_abundance_plot.py \
    -i /Metrico_project/input_file/all_cog_matrix.csv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output --prefix COG
```

### D. EC Matrix Heatmap
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/ec_heatmap.py \
    -i /Metrico_project/input_file/all_ec_matrix.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -n 50 -o /Metrico_project/output --prefix EC
```

### E. EC Circos Plot
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/ec_circos_plot.py \
    -i /Metrico_project/input_file/all_sample_ec.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output --prefix EC_Circos
```

### F. KO Matrix Heatmap
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/ko_heatmap.py \
    -i /Metrico_project/input_file/all_ko_matrix.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -n 50 -o /Metrico_project/output --prefix KO
```

### G. KO Circos Plot
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/ko_circos_plot.py \
    -i /Metrico_project/input_file/all_sample_ko.tsv \
    -mod /Metrico_project/input_file/all_module.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output --prefix KO_Circos
```

### H. MetaCyc Pathway Bubble Plot
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/pathway_bubble_plot.py \
    -i /Metrico_project/input_file/all-pathabundance_cpm.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output --prefix MetaCyc
```

### I. Ordination & PERMANOVA Analysis
```bash
sudo docker run --rm -v $(pwd):/Metrico_project mantrilabnabi/metrico_viz:v1 \
  python3 /Metrico_project/ordination_permanova.py \
    -i /Metrico_project/input_file/all_ko_matrix.tsv \
    -m /Metrico_project/input_file/Plant-host.csv \
    -o /Metrico_project/output --prefix Ordination \
    --metric aitchison
```

---

## 📊 Output File Deliverables

Each script outputs high-resolution static PNG figures suitable for publication alongside interactive HTML files for browser exploration:

* **CAZy / GH / EC / KO Heatmaps:** Generates Z-score or relative abundance heatmaps (`<PREFIX>.png` and `<PREFIX>.html`).
* **COG Module:** Generates individual functional bar plots (`<PREFIX>_<CATEGORY>_relative_abundance.png` and `.html`) plus a summary heatmap (`<PREFIX>_ALL_COG_categories_summary_heatmap.html`).
* **EC & KO Circos Modules:** Produce circular genomic/metabolic chord plots (`<PREFIX>_plot.png`) and interactive shared feature matrices (`<PREFIX>_shared_matrix.html`).
* **Pathway Bubble Module:** Produces faceted bubble plots mapping pathway abundance and peak taxon drivers (`<PREFIX>.png` and `<PREFIX>.html`).
* **Ordination & PERMANOVA Module:** Produces PCoA scatter plots with 95% confidence ellipses (`<PREFIX>_pcoa_plot.png`) and summary statistics tables (`<PREFIX>_permanova_results.txt`).

---

## 📧 Contact & Support

For queries, bug reports, feature requests, or collaboration regarding Metrico_project, please contact:

    Shrikant Mantri: shrikant@nabi.res.in

    Ardhendu Chakrabortty: ardhenduchakraborty18@gmail.com

    Computational Biology Laboratory (CBL)

    National Agri-Food Biotechnology Institute (NABI), Mohali, Punjab, India
