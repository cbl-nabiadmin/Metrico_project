#!/usr/bin/env python3
"""
Master Pipeline Runner for Metrico Visualization Scripts
--------------------------------------------------------
Executes all functional visualization scripts (CAZy, GH, COG, EC, KO, 
Pathway Bubble Plots, and Ordination/PERMANOVA) sequentially.
"""

import os
import sys
import argparse
import subprocess
import time

# =========================================================
# 1. COMMAND LINE ARGUMENTS
# =========================================================
parser = argparse.ArgumentParser(
    description="Master Pipeline Runner for Metrico Visualization Scripts",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Global settings across all pipeline scripts
parser.add_argument(
    "-m", "--metadata", 
    type=str, 
    default="Plant-host.csv",
    help="Filename or path of the metadata CSV file (inside input_file/)"
)
parser.add_argument(
    "-o", "--output-dir", 
    type=str, 
    default="/Metrico_project/output",
    help="Path to the output directory where plots and stats will be saved"
)
parser.add_argument(
    "--prefix", 
    type=str, 
    default="Metrico",
    help="Global prefix appended to generated output files"
)

# Module File Input (Explicit dest="mod" guarantees args.mod availability)
parser.add_argument(
    "--mod", "--input-module", 
    dest="mod",
    type=str, 
    default="all_module.tsv",
    help="Filename or relative path of the Module long-format file (inside input_file/)"
)

# Script-specific input overrides
parser.add_argument(
    "--cazy-class-input", 
    type=str, 
    default="cazy_class_matrix.tsv",
    help="Input file for CAZy Class heatmap"
)
parser.add_argument(
    "--gh-family-input", 
    type=str, 
    default="gh_family_filtered.csv",
    help="Input file for GH Family heatmap"
)
parser.add_argument(
    "--cog-input", 
    type=str, 
    default="all_cog_matrix.csv",
    help="Input file for COG barplots"
)
parser.add_argument(
    "--ec-matrix-input", 
    type=str, 
    default="all_ec_matrix.tsv",
    help="Input file for EC Matrix heatmap"
)
parser.add_argument(
    "--ec-circos-input", 
    type=str, 
    default="all_sample_ec.tsv",
    help="Input file for EC Circos plot"
)
parser.add_argument(
    "--ko-matrix-input", 
    type=str, 
    default="all_ko_matrix.tsv",
    help="Input file for KO Matrix heatmap"
)
parser.add_argument(
    "--ko-circos-input", 
    type=str, 
    default="all_sample_ko.tsv",
    help="Input file for KO Circos plot"
)
parser.add_argument(
    "--pathway-input", 
    type=str, 
    default="all-pathabundance_cpm.tsv",
    help="Input file for Pathway bubble plot"
)

# Ordination & PERMANOVA specific inputs
parser.add_argument(
    "--ordination-input", 
    type=str, 
    default="all_ko_matrix.tsv",
    help="Input wide matrix for PCoA & PERMANOVA (defaults to --ko-matrix-input if unset)"
)
parser.add_argument(
    "--ordination-metric", 
    type=str, 
    choices=["aitchison", "bray"], 
    default="aitchison",
    help="Distance metric for PCoA ordination (aitchison or bray)"
)
parser.add_argument(
    "--ellipse-min-n", 
    type=int, 
    default=3,
    help="Minimum group size to compute and draw 95%% group confidence ellipses"
)
parser.add_argument(
    "--permutations", 
    type=int, 
    default=999,
    help="Number of permutations for PERMANOVA hypothesis testing"
)

args = parser.parse_args()

# =========================================================
# 2. RESOLVE DIRECTORY PATHS & MODULE REGISTRATION
# =========================================================
BASE_DIR = "/Metrico_project"
OUTPUT_DIR = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(BASE_DIR, args.output_dir)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Map child scripts to their respective command-line input argument values
PIPELINE_MODULES = [
    {
        "script": "cazy_heatmap.py",
        "input": args.cazy_class_input,
        "prefix_tag": f"{args.prefix}_CAZy_Class"
    },
    {
        "script": "gh_family_heatmap.py",
        "input": args.gh_family_input,
        "prefix_tag": f"{args.prefix}_GH_Family"
    },
    {
        "script": "cog_abundance_plot.py",
        "input": args.cog_input,
        "prefix_tag": f"{args.prefix}_COG"
    },
    {
        "script": "ec_heatmap.py",
        "input": args.ec_matrix_input,
        "prefix_tag": f"{args.prefix}_EC_Matrix"
    },
    {
        "script": "ec_circos_plot.py",
        "input": args.ec_circos_input,
        "prefix_tag": f"{args.prefix}_EC_Circos"
    },
    {
        "script": "ko_heatmap.py",
        "input": args.ko_matrix_input,
        "prefix_tag": f"{args.prefix}_KO_Matrix"
    },
    {
        "script": "ko_circos_plot.py",
        "input": args.ko_circos_input,
        "module_input": args.mod,  # Forwards module mapping file to KO circos
        "prefix_tag": f"{args.prefix}_KO_Circos"
    },
    {
        "script": "pathway_bubble_plot.py",
        "input": args.pathway_input,
        "prefix_tag": f"{args.prefix}_Pathway_Bubble"
    },
    {
        "script": "ordination_permanova.py",
        "input": args.ordination_input if args.ordination_input else args.ko_matrix_input,
        "extra_flags": [
            "--metric", args.ordination_metric,
            "--ellipse-min-n", str(args.ellipse_min_n),
            "--permutations", str(args.permutations)
        ],
        "prefix_tag": f"{args.prefix}_Ordination"
    }
]

def run_script(module):
    """Runs a single pipeline script passing the specified CLI arguments."""
    script_name = module["script"]
    script_path = os.path.join(BASE_DIR, script_name)
    
    print("\n" + "=" * 60)
    print(f"🚀 RUNNING: {script_name}")
    print("=" * 60)
    
    if not os.path.exists(script_path):
        print(f"⚠️ Warning: Script '{script_path}' not found. Skipping...")
        return False

    # Base execution command
    cmd = [
        sys.executable, script_path,
        "-i", module["input"],
        "-m", args.metadata,
        "-o", OUTPUT_DIR,
        "--prefix", module["prefix_tag"]
    ]

    # Explicitly append --input-module if child script accepts it
    if "module_input" in module:
        cmd.extend(["--input-module", module["module_input"]])
        print(f"📂 Forwarded Module File:  {module['module_input']}")

    # Append module-specific extra flags if present (e.g. for ordination)
    if "extra_flags" in module:
        cmd.extend(module["extra_flags"])
        print(f"⚙️ Forwarded Extra Flags:  {' '.join(module['extra_flags'])}")

    print(f"📂 Forwarded Input File:   {module['input']}")
    print(f"📂 Forwarded Metadata:     {args.metadata}")
    print(f"📂 Output Directory:       {OUTPUT_DIR}")
    print(f"🏷️  Output Prefix Tag:     {module['prefix_tag']}")

    start_time = time.time()
    
    # Execute child script via subprocess
    result = subprocess.run(cmd, capture_output=False)
    
    elapsed = time.time() - start_time
    
    if result.returncode == 0:
        print(f"✅ PASSED: {script_name} ({elapsed:.2f} seconds)")
        return True
    else:
        print(f"❌ FAILED: {script_name} (Exit code: {result.returncode})")
        return False

def main():
    print("=" * 60)
    print("🌟 STARTING METRICO VISUALIZATION PIPELINE")
    print("=" * 60)
    print(f"📁 Working Directory: {BASE_DIR}")
    print(f"📁 Output Directory:  {OUTPUT_DIR}")
    print(f"🏷️  Global Prefix:     {args.prefix}")
    print(f"📊 Total Modules:      {len(PIPELINE_MODULES)}")

    pipeline_start = time.time()
    passed = 0
    failed = 0
    failed_scripts = []

    for module in PIPELINE_MODULES:
        success = run_script(module)
        if success:
            passed += 1
        else:
            failed += 1
            failed_scripts.append(module["script"])

    total_time = time.time() - pipeline_start
    
    print("\n" + "=" * 60)
    print("📊 PIPELINE SUMMARY REPORT")
    print("=" * 60)
    print(f"⏱️ Total Execution Time: {total_time / 60:.2f} minutes")
    print(f"✅ Successful: {passed}/{len(PIPELINE_MODULES)}")
    print(f"❌ Failed:     {failed}/{len(PIPELINE_MODULES)}")
    
    if failed_scripts:
        print("\nFailed scripts:")
        for fs in failed_scripts:
            print(f"  - {fs}")
        sys.exit(1)
    else:
        print("\n🎉 ALL VISUALIZATIONS GENERATED SUCCESSFULLY!")
        print(f"📂 Output files saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
