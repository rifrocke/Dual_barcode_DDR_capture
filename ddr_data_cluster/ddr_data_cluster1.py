import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
import warnings

# Suppress FutureWarning from pandas.concat (though it's good practice to address them specifically)
warnings.simplefilter(action='ignore', category=FutureWarning)

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description='Process DDR data from an Excel file and generate statistics and a heatmap.')

parser.add_argument("-i_excel", dest="input_excel_path", type=str, required=True,
                    help="Path to the input Excel file (e.g., analyzed_U6DDR_naive_1_w_annotation.xlsx).")
parser.add_argument("-i_csv_heatmap", dest="input_csv_heatmap_path", type=str,
                    help="Optional: Path to the input CSV file for the heatmap (e.g., U6DDR_Cycling_CLV_pct_average.csv). If not provided, the script will use the 'stats_pct' generated from the Excel file for the heatmap.")
parser.add_argument("-o", dest="output_dir", type=str, default="DDR_analysis_output",
                    help="Output directory for CSV stats and heatmap image (Default: 'DDR_analysis_output' in current location).")

args = parser.parse_args()

# --- Setup Output Directory ---
# Create the output directory if it doesn't exist.
os.makedirs(args.output_dir, exist_ok=True)
print(f"Output directory set to: {os.path.abspath(args.output_dir)}")

# --- Load Data from Excel ---
try:
    df = pd.read_excel(args.input_excel_path)
    print(f"Successfully loaded data from: {args.input_excel_path}")
    print("First 5 rows of the loaded DataFrame:")
    print(df.head(5))
except FileNotFoundError:
    print(f"❌ Error: Input Excel file not found at {args.input_excel_path}")
    exit()
except Exception as e:
    print(f"❌ Error loading Excel file: {e}")
    exit()

# Select relevant columns for analysis
df_s = df[['ID', 'Symbol', 'D_length', 'I_length']].copy() # Use .copy() to avoid SettingWithCopyWarning
print("\nSelected columns for analysis:")
df_s.info()

# Basic data checks
unique_symbols_count = len(df_s["Symbol"].unique())
print(f"\nNumber of unique 'Symbol' values: {unique_symbols_count}")

null_lengths_count = len(df_s.loc[(df_s['D_length'].isnull() & df_s['I_length'].isnull())])
print(f"Number of rows with both 'D_length' and 'I_length' as null: {null_lengths_count}")

# --- Calculate Statistics ---
stats = pd.DataFrame([])
stats_pct = pd.DataFrame([])

print("\nCalculating statistics for each gene...")
for i, gene_symbol in enumerate(df_s["Symbol"].unique()):
    gene = str(gene_symbol)
    df_sub = df_s.loc[(df_s['Symbol'] == gene)]

    # Count occurrences for different categories
    # Ensure D_length is numeric before comparison, coercing errors to NaN
    df_sub['D_length_numeric'] = pd.to_numeric(df_sub['D_length'], errors='coerce')

    Micro_D = len(df_sub.loc[(df_sub['D_length_numeric'] < 3)])
    Mid_D = len(df_sub.loc[(df_sub['D_length_numeric'] >= 3) & (df_sub['D_length_numeric'] <= 10)])
    Large_D = len(df_sub.loc[(df_sub['D_length_numeric'] > 10)])
    Ins = df_sub["I_length"].notna().sum() # Count non-null insertions
    Direct = len(df_sub.loc[(df_sub['D_length'].isnull() & df_sub['I_length'].isnull())]) # Original logic for Direct

    # Create row for absolute counts
    current_stats_row = pd.DataFrame({
        "gene": [gene],
        "Large_D(>10bps)": [Large_D],
        "Mid_D(3-10bps)": [Mid_D],
        "Micro_D(<3bps)": [Micro_D],
        "Direct(or_noDSB)": [Direct],
        "Insertion": [Ins]
    })
    stats = pd.concat([stats, current_stats_row], ignore_index=True)

    # Calculate percentages
    Total = Micro_D + Mid_D + Large_D + Ins + Direct
    if Total == 0: # Avoid division by zero
        Micro_D_pct = Mid_D_pct = Large_D_pct = Ins_pct = Direct_pct = 0.0
    else:
        # Note: Original code had negative percentages for deletions. Keeping that logic.
        Micro_D_pct = -Micro_D / Total
        Mid_D_pct = -Mid_D / Total
        Large_D_pct = -Large_D / Total
        Ins_pct = Ins / Total
        Direct_pct = Direct / Total

    # Create row for percentage counts
    current_stats_pct_row = pd.DataFrame({
        "gene": [gene],
        "Large_D(>10bps)": [Large_D_pct],
        "Mid_D(3-10bps)": [Mid_D_pct],
        "Micro_D(<3bps)": [Micro_D_pct],
        "Direct(or_noDSB)": [Direct_pct],
        "Insertion": [Ins_pct]
    })
    stats_pct = pd.concat([stats_pct, current_stats_pct_row], ignore_index=True)

print("Statistics calculation complete.")

# --- Save Statistics to CSV ---
output_csv_path_stats = os.path.join(args.output_dir, 'U6DDR_naive_1.csv')
output_csv_path_stats_pct = os.path.join(args.output_dir, 'U6DDR_naive_1_pct.csv')

stats.to_csv(output_csv_path_stats, index=False)
stats_pct.to_csv(output_csv_path_stats_pct, index=False)
print(f"Absolute counts saved to: {output_csv_path_stats}")
print(f"Percentage counts saved to: {output_csv_path_stats_pct}")

# --- Prepare Data for Heatmap ---
# If a specific CSV for heatmap is provided, use it. Otherwise, use the 'stats_pct' generated above.
if args.input_csv_heatmap_path:
    try:
        stats_pct_for_heatmap = pd.read_csv(args.input_csv_heatmap_path)
        print(f"\nUsing provided CSV for heatmap: {args.input_csv_heatmap_path}")
    except FileNotFoundError:
        print(f"❌ Error: Heatmap input CSV file not found at {args.input_csv_heatmap_path}. Using generated stats_pct for heatmap.")
        stats_pct_for_heatmap = stats_pct.copy()
    except Exception as e:
        print(f"❌ Error loading heatmap CSV: {e}. Using generated stats_pct for heatmap.")
        stats_pct_for_heatmap = stats_pct.copy()
else:
    print("\nNo specific CSV for heatmap provided. Using the 'stats_pct' generated from the Excel file.")
    stats_pct_for_heatmap = stats_pct.copy()

# Set 'gene' as index for the heatmap
if 'gene' in stats_pct_for_heatmap.columns:
    stats_pct_for_heatmap.set_index("gene", inplace=True)
else:
    print("⚠️ Warning: 'gene' column not found in heatmap data. Heatmap might not have gene labels.")

print("\nHeatmap data info:")
stats_pct_for_heatmap.info()
print("First 5 rows of heatmap data:")
print(stats_pct_for_heatmap.head(5))

# --- Generate and Save Heatmap ---
print("\nGenerating heatmap...")
# Adjust figure size for better visibility, especially with many genes
plt.figure(figsize=(10, max(6, len(stats_pct_for_heatmap) * 0.15))) # Dynamic height based on number of genes

sns_plot = sns.clustermap(stats_pct_for_heatmap,
                          cmap="coolwarm",
                          yticklabels=True,
                          col_cluster=False, # Columns are not clustered
                          cbar_pos=(0.02, 0.85, 0.02, 0.18), # Position of the color bar
                          figsize=(10, max(6, len(stats_pct_for_heatmap) * 0.15)) # Ensure figsize is also passed here
                         )

# Adjust heatmap and dendrogram positions for better layout
# This part of the code is highly specific to the original notebook's layout
# and might need manual tweaking depending on the number of genes and desired output.
hmap = sns_plot.ax_heatmap.get_position()
# Adjust heatmap width (e.g., 25% of original width)
sns_plot.ax_heatmap.set_position([hmap.x0, hmap.y0, hmap.width * 0.25, hmap.height])

# Adjust column dendrogram position (if col_cluster was True, otherwise it's just a placeholder)
col = sns_plot.ax_col_dendrogram.get_position()
sns_plot.ax_col_dendrogram.set_position([col.x0, col.y0, col.width * 0.25, col.height * 0.5])

# Adjust y-axis tick label font size
plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)

# Save the heatmap image
output_png_path = os.path.join(args.output_dir, "U6DDR_Cycling_V2.png")
sns_plot.savefig(output_png_path, dpi=600, bbox_inches='tight') # bbox_inches='tight' helps prevent labels from being cut off
print(f"Heatmap saved to: {output_png_path}")

# Display the plot (optional for local execution, useful for interactive sessions)
plt.show()