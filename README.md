🧬Metrico_project : Metagenomic Functional Visualization Pipeline
Welcome to Metrico_project! This toolkit provides a streamlined Python/Docker pipeline for processing and visualizing metagenomic functional annotations across plant hosts and environmental samples (including CAZy Classes, GH Families, COG Categories, EC Matrices, KO Matrices, EC Circos Plots, KO Circos Plots, and MetaCyc Pathway Bubble Plots). The Docker container comes equipped with an input_file/ directory containing example files to help you quickly structure your data, alongside pre-packaged scripts for all necessary data preprocessing and visualization tasks; a complete step-by-step guide detailing how to utilize these scripts is provided below in this document.
📌 Navigation Index
🧭 Workflow Fast-Track & Intermediate Entry Points
⚙️ Prerequisites & Data Pre-Processing Tutorial
📁 Project Directory Structure
🛠️ Step 1: Unpack & Load the Docker Container
🚀 Step 2: Running the Master Pipeline (run_all.py)
📋 Full Parameter Reference Table for run_all.py
🎯 Running Individual Scripts Separately
📊 Output File Deliverables
🧭 Workflow Fast-Track & Intermediate Entry Points
You do not need to run every pre-processing step or tool if you already have processed matrices or only want specific plots. Use the table below to determine where to start:
Your Current Starting Point
Can You Skip Pre-Processing?
Actions Required
I want to test the pipeline using provided example data
YES (Skip Tasks 1–5)
Jump directly to Step 1 and run run_all.py with default settings.
I already have formatted matrices (e.g., all_ko_matrix.tsv, cazy_class_matrix.tsv)
YES (Skip Tasks 1–5)
Place your files into the input_file/ directory and jump directly to Step 2 or Individual Scripts.
I only want to run specific plots (e.g., only CAZy or COG)
PARTIAL
Only execute the specific pre-processing tasks required for your target plot then run the corresponding Individual Script.
I have raw eggNOG (.emapper.annotations) or HUMAnN (pathabundance.tsv) outputs
NO
Complete the relevant pre-processing tasks below before running any visualization scripts.

⚙️ Prerequisites & Data Pre-Processing Tutorial
⚠️ CRITICAL REQUIREMENT: Unless you already have pre-formatted matrices, you must first process your raw FASTQ files using HUMAnN and your assembled FASTA files using eggNOG-mapper (following assembly with SPAdes). Once processed, convert the resulting outputs into the required matrix formats. Refer to the tutorial below to generate compatible input files inside your input_file/ directory.
Pre-Processing Workflow Map
Task
Input File
Target Visualization
Required Output Name
Combine Annotations
*.emapper.annotations
Central Master Dataset
all_with_source.annotations
Circos KO & Modules
all_with_source.annotations
Circular KO/Module Plot
all_sample_ko.tsv, all_module.tsv
Circos EC Numbers
all_with_source.annotations
Circular EC Chord Diagram
all_sample_ec.tsv
HUMAnN Abundance
*_pathabundance.tsv
Pathway Driver Bubble Plot
all-pathabundance_cpm.tsv
COG Category Matrix
all_with_source.annotations
COG Functional Barplot
all_cog_matrix.csv
KO Matrix
all_with_source.annotations
Functional KO Heatmap
all_ko_matrix.tsv
EC Matrix
all_with_source.annotations
Enzyme EC Heatmap
all_ec_matrix.tsv
CAZy Profile Matrix
all_with_source.annotations
CAZy Class/Family Heatmap
cazy_class_matrix.tsv,
gh_family_filtered.csv / gt_family_filtered.csv

Step-by-Step Pre-Processing Commands
Task 1: Combine eggNOG Annotations
Merge individual sample subfolder .emapper.annotations files into a single master annotation file with sample labels attached. Save this script as combined_results.sh and run it within the directory where emapper output files for each sample are created.
#!/bin/bash
OUT="all_with_source.annotations"
DIR="/Path/to/emapper-result/samplewise"

# Get first annotation file (from any subdirectory) to extract header
first_file=$(find "$DIR" -type f -name "*.emapper.annotations" | head -n 1)

# Extract header and add Sample column
echo -e "Sample\t$(grep "^#query" "$first_file" | sed 's/^#//')" > "$OUT"

# Loop through each subdirectory
for d in "$DIR"/*/; do
    sample=$(basename "$d")

    # Find annotation file inside that directory
    f=$(find "$d" -maxdepth 1 -name "*.emapper.annotations")

    # Skip if no file found
    [ -z "$f" ] && continue

    # Append data with sample name
    grep -v "^#" "$f" | awk -v s="$sample" '{print s"\t"$0}'
done >> "$OUT"


Run:
bash combined_results.sh
Output: all_with_source.annotations

Task 2: Circos Inputs (KO & Modules, EC Numbers)
#Extract KO IDs and unroll comma-separated KEGG Modules into long-format files. Extract sample IDs and annotated Enzyme Commission (EC) numbers. Save this script as circos_script.sh and run it within the directory where all_with_source.annotations file was created.
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

#Extract sample IDs and annotated Enzyme Commission (EC) numbers.
awk -F '\t' '$12 != "" && $12 != "-" {print $1"\t"$12}' all_with_source.annotations > all_sample_ec.tsv
Run:
bash circos_script.sh
Output: all_sample_ko.tsv, all_module.tsv, all_sample_ec.tsv

Task 3: HUMAnN Pathway CPM Normalization
Merge per-sample HUMAnN pathway abundances and normalize to Counts Per Million (CPM). Save this script as humann_rejoin.sh and run it within the directory where HUMAnN output files are created.
#!/bin/bash
# Activate HUMAnN environment before running
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
# Cleanup staging directory after processing
rm -rf "$JOIN_DIR"
Run:
bash humann_rejoin.sh
Output: all-pathabundance_cpm.tsv
Task 4: Build COG, EC & KO Matrix
#Extract primary COG category letters and shape into a Sample x COG matrix for bar plots. Extract KO IDs, expand multi-assigned entries, and format into a KO x Sample matrix. Parse EC numbers and format into an EC x Sample matrix for heatmaps. Extract CAZy annotations, derive parent CAZy classes (GH, GT, PL, etc.), and build long-format count tables.
#!/bin/bash
#Extract primary COG category letters and shape into a Sample x COG matrix for bar plots.
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

#Extract KO IDs, expand multi-assigned entries, and format into a KO x Sample matrix.
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

#Parse EC numbers and format into an EC x Sample matrix for heatmaps.
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

#Extract CAZy annotations, derive parent CAZy classes (GH, GT, PL, etc.), and build long-format count tables.
awk -F '\t' 'NR>1 { if($20 != "-" && $20 != "") print $1 "\t" $20 }' all_with_source.annotations | \
awk -F '\t' '{
    n=split($2, a, ",");
    for(i=1;i<=n;i++){ gsub(/^ +| +$/, "", a[i]); if(a[i] != "" && a[i] != "-") print $1 "\t" a[i]; }
}' > all_cazy_expanded.tsv

awk -F '\t' '{ class=$2; sub(/[0-9]+.*/, "", class); print $1 "\t" $2 "\t" class; }' all_cazy_expanded.tsv > all_cazy_class.tsv
export LC_ALL=C
sort all_cazy_class.tsv | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > all_cazy_count.csv

Run:
bash circos_script.sh
Output: all_cog_matrix.csv, all_ko_matrix.tsv, all_ec_matrix.tsv, all_cazy_count.csv 

Task 5: CAZy Profile and CAZy subfamilies Matrix
#convert all_cazy_count table to the desired matrix format using a python script, save this script as cazy_matrix_convert.py
#!/usr/bin/env python3
import pandas as pd
import sys
# File paths
input_file = "all_cazy_count.csv"  # Your long-format input file
output_file = "cazy_class_matrix.tsv"     # Generated wide-matrix file

# 1. Read input long-format table
try:
    # Handles headerless or headered 4-column files cleanly
    df = pd.read_csv(input_file, sep=None, engine="python", header=None)
    if df.iloc[0, 0].lower().startswith("sample"):
        df = df.iloc[1:]  # Drop header row if present
    df.columns = ["SampleID", "Family", "Class", "Count"]
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)

# 2. Clean data types
df["SampleID"] = df["SampleID"].astype(str).str.strip()
df["Class"] = df["Class"].astype(str).str.strip()
df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0)

# 3. Aggregate counts by (SampleID, Class) and pivot to wide matrix
matrix_df = df.groupby(["SampleID", "Class"])["Count"].sum().unstack(fill_value=0)

# 4. Save as TSV file with 'Sample' as the index header name
matrix_df.index.name = "Sample"
matrix_df.to_csv(output_file, sep="\t")

print(f"✅ Conversion complete! Matrix saved to: {output_file}")
print(matrix_df.head())
Run:
python3 cazy_matrix_convert.py
Output: cazy_class_matrix.tsv
#filter CAZy classes of your preference (GH, GT, CBM etc.), save this script as - filter_cazy_class.py
#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Generic CAZy Class Extractor",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-i", "--input", type=str, required=True, help="Input raw CAZy CSV/TSV file")
parser.add_argument("-c", "--cazy-class", type=str, default="GH", help="Target CAZy class to extract (e.g., GH, AA, CBM, CE, PL, GT)")
parser.add_argument("-o", "--output", type=str, default=None, help="Output CSV path (optional)")

args = parser.parse_args()

target_class = args.cazy_class.upper().strip()
output_file = args.output if args.output else f"{target_class.lower()}_family_filtered.csv"

print(f"📂 Loading input file: {args.input}")
print(f"🎯 Target CAZy Class: {target_class}")

# 1. Read input data handling headered/headerless formats
try:
    df = pd.read_csv(args.input, sep=None, engine="python", header=None)
    if str(df.iloc[0, 0]).lower().startswith("sample"):
        df = df.iloc[1:]  # Drop header row if present
    df = df.iloc[:, :4]
    df.columns = ["SampleID", "Family", "Class", "Count"]
except Exception as e:
    print(f"❌ Error loading file: {e}")
    sys.exit(1)

# 2. Clean data fields
df["SampleID"] = df["SampleID"].astype(str).str.strip()
df["Family"] = df["Family"].astype(str).str.strip()
df["Class"] = df["Class"].astype(str).str.upper().str.strip()
df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0).astype(int)

# 3. Filter dynamically for the requested class or family prefix
filtered_df = df[(df["Class"] == target_class) | (df["Family"].str.startswith(target_class))].copy()

if filtered_df.empty:
    print(f"⚠️ Warning: No entries found for class '{target_class}'!")
    sys.exit(0)

# 4. Format to target column arrangement: [Count, SampleID, <CLASS>_Family]
family_col_name = f"{target_class}_Family"
filtered_df = filtered_df.rename(columns={"Family": family_col_name})
filtered_df = filtered_df[["Count", "SampleID", family_col_name]]

# 5. Export result
filtered_df.to_csv(output_file, index=False)
print(f"✅ Extracted {len(filtered_df)} entries for '{target_class}'. Saved to: {output_file}")
print(filtered_df.head(10))
Run:  
python3 filter_cazy_class.py -i all_cazy_count.csv -c GH <class name (GH/GT/CBM etc.)> -o gh_family_filtered.csv
Output: gh_family_filtered.csv

📁 Project Directory Structure
Once pre-processing is complete, place your generated files into the input_file/ directory as structured below:

Metrico_project/
├── input_file/                       # Store input matrices & metadata here
│   ├── cazy_class_matrix.tsv         # Example CAZy Class matrix
│   ├── gh_family_filtered.csv        # Example GH Family count matrix
│   ├── all_cog_matrix.csv           # Example COG functional matrix
│   ├── all_ec_matrix.tsv             # Example EC count matrix
│   ├── all_sample_ec.tsv          # Example long-format EC file
│   ├── all_ko_matrix.tsv             # Example KO count matrix
│   ├── all_sample_ko.tsv          # Example long-format KO file
│   ├── all-pathabundance_cpm.tsv     # Example MetaCyc pathway abundance file
│   ├── MetaCyc-pathways.txt          # MetaCyc ontology mapping file taken from MetaCyc database portal
│   └── Plant-host.csv               # Example plant host metadata
├── output/                           # Results directory
├── cazy_heatmap.py             # CAZy Class analysis script
├── gh_family_heatmap.py              # GH Family analysis script
├── cog_abundance_plot.py     # COG analysis script
├── ec_heatmap.py              # EC Matrix analysis script
├── ec_circos_plot.py                 # EC Circos analysis script
├── ko_heatmap.py                     # KO Matrix analysis script
├── ko_circos_plot.py                 # KO Circos analysis script
├── pathway_bubble_plot.py            # Pathway Bubble plot script
└── run_all.py                        # Master pipeline wrapper script
🛠️ Step 1: Unpack & Load the Docker Container
The pipeline environment is delivered as a compressed archive (.tar.gz). Before running any script, load the image into Docker:

Navigate to the folder containing your .tar.gz file:
cd /path/to/image/archive
Load the Docker image:
sudo docker load -i metrico_viz:v1.tar.gz
Verify that the image has loaded successfully:
sudo docker images
(You should see metrico_viz:v1 listed under REPOSITORY).
🚀 Step 2: Running the Master Pipeline (run_all.py)
1. Default Run (Example Data)
Example input datasets are pre-packaged inside the input_file/ directory. Running run_all.py without extra flags processes all default example datasets sequentially and outputs figures into output/.
Bash
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/run_all.py

2. Full Master Command (Custom Data Across All 9 Scripts)
To run custom input datasets across all 9 modules without editing any Python scripts, pass custom file names directly via command-line arguments:
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/run_all.py \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix> \
    --mod <path to all_module.tsv file> \
    --cazy-class-input <path to CAZy class matrix file> \
    --gh-family-input <path to GH family count file> \
    --cog-input <path to COG functional matrix file> \
    --ec-matrix-input <path to EC matrix file> \
    --ec-circos-input <path to long-format EC file> \
    --ko-matrix-input <path to KO matrix file> \
    --ko-circos-input <path to long-format KO file> \
    --pathway-input <path to pathway abundance file> \
    --ordination-input <path to KO matrix file> \
    --ordination-metric <aitchison / bray> \
    --ellipse-min-n 3 \
    --permutations 999
    
📋 Full Parameter Reference Table for run_all.py
Parameter Flag
Description
Default Value (inside input_file/)
--rm
Container Cleanup Flag: Automatically removes the anonymous Docker container filesystem upon exit to prevent disk clutter.
Docker Runtime Flag
-m / --metadata
Common metadata CSV file, for example Plant-host.csv file is provided
Plant-host.csv
-o / --output-dir
Target folder location for all generated output files
/Metrico_project/output
--prefix
Global prefix tag prepended to saved output filenames
Metrico
--mod    / --input-module
Module file required for KO circos plot
all_module.tsv
--cazy-class-input
Input file for cazy_heatmap.py
cazy_class_matrix.tsv
--gh-family-input
Input file for gh_family_heatmap.py
gh_family_filtered.csv / gt_family_filtered.csv
--cog-input
Input file for cog_abundance_plot.py
all_cog_matrix.csv
--ec-matrix-input
Input file for ec_heatmap.py
all_ec_matrix.tsv
--ec-circos-input
Input file for ec_circos_plot.py
all_sample_ec.tsv
--ko-matrix-input
Input file for ko_heatmap.py
all_ko_matrix.tsv
--ko-circos-input
Input file for ko_circos_plot.py
all_sample_ko.tsv
--pathway-input
Input file for pathway_bubble_plot.py
all-pathabundance_cpm.tsv
--ordination-input
Input file for ordination_permanova.py 
all_ko_matrix.tsv
--ordination-metric
Distance metric for ordination (`aitchison` = CLR + Euclidean, `bray` = Bray-Curtis)
aitchison
--ellipse-min-n
Minimum host group size required to compute and draw 95% confidence ellipses
3
--permutations
Number of permutations for PERMANOVA hypothesis testing
999

🎯 Running Individual Scripts Separately
You can also execute any analysis script independently using -i (input), -m (metadata), -o (output directory), and --prefix (output file prefix):
A. CAZy Class Heatmap 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/cazy_heatmap.py \
    -i <path to CAZy class matrix file> \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix>

B. GH Family Heatmap 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/gh_family_heatmap.py \
    -i <path to GH family count file> \
    -m <path to metadata file> \
    -n <no. of gh family to be plotted> \
    -o <path to output directory> \
    --prefix <output file prefix>

C. COG Categories Bar Plots & Heatmap 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/cog_abundance_plot.py \
    -i <path to COG functional matrix file> \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix>

D. EC Matrix Heatmap 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/ec_heatmap.py \
    -i <path to EC matrix file> \
    -m <path to metadata file> \
    -n <no. of ECs to be plotted> \
    -o <path to output directory> \
    --prefix <output file prefix>

E. EC Circos Plot 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/ec_circos_plot.py \
    -i <path to long-format EC file> \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix>

F. KO Matrix Heatmap 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/ko_heatmap.py \
    -i <path to KO matrix file> \
    -m <path to metadata file> \
    -n <no. of KOs to be plotted> \
    -o <path to output directory> \
    --prefix <output file prefix>

G. KO Circos Plot 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/ko_circos_plot.py \
    -i <path to long-format KO file> \
    --mod <path to module file> \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix>

H. MetaCyc Pathway Bubble Plot 
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/pathway_bubble_plot.py \
    -i <path to pathway abundance file> \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix>

I. Ordination & PERMANOVA Analysis
sudo docker run --rm \
  -v /path/to/Metrico_project:/Metrico_project \
  metrico_viz:v1 python3 /Metrico_project/ordination_permanova.py \
    -i <path to KO matrix file> \
    -m <path to metadata file> \
    -o <path to output directory> \
    --prefix <output file prefix> \
    --metric < aitchison / bray> \


📊 Output File Deliverables
Each script outputs high-resolution static PNG figures suitable for publication alongside interactive HTML files for browser exploration:
CAZy / GH / EC / KO Heatmaps: Generates Z-score or relative abundance heatmaps (<PREFIX>.png and <PREFIX>.html).

COG Module: Generates individual functional bar plots (<PREFIX>_<CATEGORY>_relative_abundance.png and .html) plus a summary heatmap (<PREFIX>_ALL_COG_categories_summary_heatmap.html).

EC & KO Circos Modules: Produce circular genomic/metabolic chord plots (<PREFIX>_plot.png) and interactive shared feature matrices (<PREFIX>_shared_matrix.html).



Pathway Bubble Module: Produces faceted bubble plots mapping pathway abundance and peak taxon drivers (<PREFIX>.png and <PREFIX>.html).








