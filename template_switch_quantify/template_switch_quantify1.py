import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os # Import os for path manipulation
import argparse # Import argparse for command-line arguments

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description='Analyze dual-barcode library data, merge R1 and R2 reads, calculate alignment scores, and identify specific read types.')

parser.add_argument("-r1", dest="r1_csv_path", type=str, required=True,
                    help="Path to the R1 CSV file (e.g., p320_RuDDR_Taq3_R1.csv).")
parser.add_argument("-r2", dest="r2_csv_path", type=str, required=True,
                    help="Path to the R2 CSV file (e.g., p320_RuDDR_Taq3_R2.csv).")
parser.add_argument("-ref", dest="ref_csv_path", type=str, required=True,
                    help="Path to the reference CSV file (e.g., human_DDR_minipool_ref.csv).")
parser.add_argument("-o", dest="output_dir", type=str, default="barcode_analysis_output",
                    help="Output directory for processed CSV files (Default: 'barcode_analysis_output' in current location).")

args = parser.parse_args()

# --- Setup Output Directory ---
# Create the output directory if it doesn't exist.
os.makedirs(args.output_dir, exist_ok=True)
print(f"Output directory set to: {os.path.abspath(args.output_dir)}")

# --- Function for Reverse Complement ---
def RCP(dna):
    """
    Calculates the reverse complement of a DNA sequence.
    Handles 'N' as 'N'.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    # Reverse the sequence, get complement for each base, and join them.
    return ''.join([complement[base] for base in dna[::-1]])

# --- Load Data ---
# Use try-except blocks for robust file loading.
try:
    df_R1 = pd.read_csv(args.r1_csv_path)
    print(f"Loaded R1 data from: {args.r1_csv_path} (Rows: {len(df_R1)})")
except FileNotFoundError:
    print(f"❌ Error: R1 CSV file not found at {args.r1_csv_path}")
    exit()
except Exception as e:
    print(f"❌ Error loading R1 CSV: {e}")
    exit()

try:
    df_R2 = pd.read_csv(args.r2_csv_path)
    print(f"Loaded R2 data from: {args.r2_csv_path} (Rows: {len(df_R2)})")
except FileNotFoundError:
    print(f"❌ Error: R2 CSV file not found at {args.r2_csv_path}")
    exit()
except Exception as e:
    print(f"❌ Error loading R2 CSV: {e}")
    exit()

try:
    df_ref = pd.read_csv(args.ref_csv_path)
    print(f"Loaded reference data from: {args.ref_csv_path} (Rows: {len(df_ref)})")
    print("First 5 rows of reference data:")
    print(df_ref.head(5))
except FileNotFoundError:
    print(f"❌ Error: Reference CSV file not found at {args.ref_csv_path}")
    exit()
except Exception as e:
    print(f"❌ Error loading reference CSV: {e}")
    exit()

# --- Preprocessing R1 and R2 DataFrames ---

# Split 'SeqID' into 'ID' and 'Run' columns for both R1 and R2.
# 'expand=True' creates new columns directly.
df_R1[['ID', 'Run']] = df_R1['SeqID'].str.split(' ', expand=True)
df_R2[['ID', 'Run']] = df_R2['SeqID'].str.split(' ', expand=True)

# Drop unnecessary columns. Use .copy() to ensure operations on new DataFrame.
df_R1_processed = df_R1.drop(["SeqID", "quality", "Run"], axis=1).copy()
df_R2_processed = df_R2.drop(["SeqID", "quality", "Run"], axis=1).copy()

print("\nProcessed R1 DataFrame head:")
print(df_R1_processed.head())
print("\nProcessed R1 DataFrame info:")
df_R1_processed.info()

# Calculate sequence lengths and add as a new column 'len'.
# Convert to string first to handle potential non-string sequence entries.
df_R1_processed["len"] = df_R1_processed["sequence"].astype(str).apply(len)
df_R2_processed["len"] = df_R2_processed["sequence"].astype(str).apply(len)

# Filter reads based on length criteria.
# This removes rows where 'len' is outside the specified range.
initial_len_R1 = len(df_R1_processed)
initial_len_R2 = len(df_R2_processed)

df_R1_filtered = df_R1_processed[(df_R1_processed.len >= 17) & (df_R1_processed.len <= 19)].copy()
df_R2_filtered = df_R2_processed[(df_R2_processed.len >= 19) & (df_R2_processed.len <= 21)].copy()

print(f"\nR1 reads filtered: {initial_len_R1 - len(df_R1_filtered)} rows removed (length < 17 or > 19).")
print(f"R2 reads filtered: {initial_len_R2 - len(df_R2_filtered)} rows removed (length < 19 or > 21).")

print("\nFiltered R1 DataFrame head:")
print(df_R1_filtered.head())
print("\nFiltered R2 DataFrame head:")
print(df_R2_filtered.head())

# --- Merge R1 and R2 DataFrames ---
# Merge based on 'ID' column using an inner join to keep only common IDs.
df_merged = pd.merge(df_R1_filtered, df_R2_filtered, on=["ID"], how='inner')
print(f"\nMerged DataFrame created with {len(df_merged)} rows.")
print("Merged DataFrame head:")
print(df_merged.head(5))

# --- Filter by Reference Barcode ---
# Keep only rows where 'sequence_x' (from R1) is present in the reference 'barcode' column.
df_clean = df_merged[df_merged['sequence_x'].isin(df_ref['barcode'])].copy()
print(f"\nCleaned DataFrame (R1 barcode in reference) has {len(df_clean)} rows.")

# --- Calculate Reverse Complement of R2 Sequence ---
# Apply the RCP function to 'sequence_y' (from R2) to get its reverse complement.
# Ensure 'sequence_y' is treated as string to avoid errors with non-string types.
print("\nCalculating reverse complement for R2 sequences...")
df_clean["sequence_y_rcp"] = df_clean["sequence_y"].astype(str).apply(RCP)

# Reset index after filtering to ensure a continuous index for iteration.
df_clean.reset_index(drop=True, inplace=True) # Use drop=True to avoid adding old index as a column

print("\nCleaned DataFrame info after RCP and index reset:")
df_clean.info()

# --- Biopython Alignment ---
# NOTE: Biopython must be installed in your local Python environment.
# If you don't have it, run: pip install biopython
try:
    from Bio import Align
    aligner = Align.PairwiseAligner()
    # Default parameters for aligner are usually good for basic sequence alignment.
    # You might want to set specific match/mismatch/gap scores if needed, e.g.:
    # aligner.match_score = 1
    # aligner.mismatch_score = -1
    # aligner.gap_score = -2
    # aligner.query_gap_score = -0.5
    # aligner.target_gap_score = -0.5
    print("\nBiopython Align module loaded.")
except ImportError:
    print("❌ Error: Biopython library not found. Please install it using 'pip install biopython'.")
    exit()
except Exception as e:
    print(f"❌ Error initializing Biopython aligner: {e}")
    exit()

print("Calculating alignment scores...")
Align_score = []
# Iterate using iterrows for safer row-wise processing, though direct column access is faster.
# For large datasets, consider vectorizing or using apply with a lambda function if possible.
for index, row in df_clean.iterrows():
    # Ensure sequences are strings before aligning
    seq_x = str(row["sequence_x"])
    seq_y_rcp = str(row["sequence_y_rcp"])
    alignments = aligner.align(seq_x, seq_y_rcp)
    # alignments[0].score gets the score of the first (best) alignment
    Align_score.append(alignments[0].score)

df_clean["Align_score"] = Align_score
print("Alignment scores calculated.")
print("Cleaned DataFrame head with Align_score:")
print(df_clean.head(5))

# --- Save Processed Data ---
output_csv_path_w_score = os.path.join(args.output_dir, 'p320_RuDDR_Taq3_w_score.csv')
df_clean.to_csv(output_csv_path_w_score, index=False)
print(f"\nProcessed data with alignment scores saved to: {output_csv_path_w_score}")

# --- Analyze Alignment Scores ---
# Count rows where 'Align_score' is greater than 16.
# Note: The original code used .count(axis=1) which counts non-NA values per row,
# but here it seems the intent was to count rows where the condition is True.
# Using sum() on a boolean series counts True values.
aligned_reads_count = (df_clean["Align_score"] > 16).sum()
print(f"\nNumber of reads with Align_score > 16: {aligned_reads_count}")

# Identify "wrong" alignments (score <= 16).
df_wrong = df_clean[df_clean["Align_score"] <= 16].copy()
print(f"Number of 'wrong' alignments (Align_score <= 16): {len(df_wrong)}")

# Identify "switch" reads: wrong alignments where RCP of R2 is in ref 'sgRNA sequence'.
# This assumes 'sgRNA sequence' is a column in your reference file.
df_switch = df_wrong[df_wrong["sequence_y_rcp"].isin(df_ref["sgRNA sequence"])].copy()
print(f"Number of 'switch' reads: {len(df_switch)}")

print("\nScript execution complete.")