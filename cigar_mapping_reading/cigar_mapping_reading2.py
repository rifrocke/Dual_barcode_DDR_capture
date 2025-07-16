import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None # Suppress SettingWithCopyWarning
import os
import argparse
import re # Import regex module for CIGAR parsing

# --- Function to Extract Mutations from CIGAR Strings ---
def mut_extraction(df_input):
    """
    Extracts deletion and insertion lengths and locations from CIGAR strings.

    Args:
        df_input (pd.DataFrame): DataFrame containing a 'Cigar' column.

    Returns:
        pd.DataFrame: A new DataFrame with 'D_length', 'D_location',
                      'I_length', and 'I_location' columns added,
                      containing lists of extracted information.
    """
    df_output = df_input.copy()
    # Initialize new columns with empty lists for each row.
    # This is crucial to ensure they are mutable lists from the start,
    # preventing TypeError when trying to append.
    df_output["D_length"] = [[] for _ in range(len(df_output))]
    df_output["D_location"] = [[] for _ in range(len(df_output))]
    df_output["I_length"] = [[] for _ in range(len(df_output))]
    df_output["I_location"] = [[] for _ in range(len(df_output))]

    # Regex to find CIGAR operations (e.g., "10M", "2I", "5D").
    # It captures the length (digits) and the operation type (letter).
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')

    # Iterate through each row of the DataFrame.
    # Using iterrows is suitable when you need to modify list-like elements
    # within DataFrame cells, though for very large datasets, vectorized
    # operations or `apply` might be more performant if the logic allows.
    for index, row in df_output.iterrows():
        cigar_str = str(row["Cigar"]) # Ensure CIGAR string is treated as a string.

        # Skip rows with no CIGAR information or NaN values.
        if cigar_str == "*" or pd.isna(cigar_str):
            continue

        # --- Custom CIGAR string parsing logic from the original code ---
        # The original code's logic (e.g., `annotation[2] == "I"`, `annotation[3:-3]`)
        # implies a very specific, non-standard CIGAR format where the actual
        # CIGAR string is embedded within a larger string with fixed prefix/suffix.
        # If your CIGAR strings are standard (e.g., "10M2I5D"), this part
        # should be removed or adjusted.
        processed_cigar_str = cigar_str
        if len(cigar_str) >= 6 and cigar_str[2] == "I":
            processed_cigar_str = cigar_str[3:-3]
        elif len(cigar_str) < 6 and cigar_str != "*":
            # Warn if string is too short for custom parsing but not empty.
            print(f"Warning: CIGAR string '{cigar_str}' at index {index} is too short for custom parsing. Treating as standard CIGAR.")

        idx0 = 0 # Initialize event location on the reference sequence for the current row.

        # Parse the CIGAR string using the regex.
        parsed_operations = cigar_pattern.findall(processed_cigar_str)

        for length_str, operation in parsed_operations:
            try:
                l = int(length_str) # Convert length string to integer.
            except ValueError:
                print(f"Warning: Could not parse length '{length_str}' in CIGAR '{cigar_str}' at index {index}. Skipping operation.")
                continue # Skip to the next operation if length is invalid.

            if operation == "M": # Match (consumes reference sequence)
                idx0 += l
            elif operation == "D": # Deletion (consumes reference sequence)
                row["D_length"].append(l)
                row["D_location"].append([idx0, idx0 + l])
                idx0 += l
            elif operation == "I": # Insertion (does NOT consume reference sequence)
                row["I_length"].append(l)
                row["I_location"].append([idx0])
            # Other CIGAR operations (N, S, H, P, =, X) are not explicitly
            # handled for length/location extraction in the original code,
            # but we need to account for their effect on `idx0` if they consume
            # reference sequence.
            elif operation in ["N", "=", "X"]: # N (skip), = (match), X (mismatch) consume reference
                idx0 += l
            # S (soft clip), H (hard clip), P (padding) do not consume reference sequence
            # for alignment, so `idx0` does not change.
            elif operation in ["S", "H", "P"]:
                pass
            else:
                print(f"Warning: Unknown CIGAR operation '{operation}' in '{cigar_str}' at index {index}. Skipping.")

    return df_output

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description='Process CIGAR strings from tab-delimited files to extract mutation information.')

parser.add_argument("-i", dest="input_folder", type=str, required=True,
                    help="Path to the input folder containing tab-delimited files with a 'Cigar' column.")
parser.add_argument("-o", dest="output_folder", type=str, default="Analyzed_U6DDR",
                    help="Path to the output folder where analyzed CSV files will be saved (Default: 'Analyzed_U6DDR').")

args = parser.parse_args()

# --- Setup Output Directory ---
# Create the output directory if it doesn't exist.
os.makedirs(args.output_folder, exist_ok=True)
print(f"Output directory set to: {os.path.abspath(args.output_folder)}")

# --- Process Files in Input Folder ---
try:
    file_names = os.listdir(args.input_folder)
    print(f"Found {len(file_names)} files in input folder: {args.input_folder}")
except FileNotFoundError:
    print(f"❌ Error: Input folder not found at {args.input_folder}")
    exit()
except Exception as e:
    print(f"❌ Error accessing input folder: {e}")
    exit()

processed_count = 0
for f in file_names:
    # Construct the full path to the input file.
    input_file_path = os.path.join(args.input_folder, f)

    # Modify this line to check for common tab-delimited extensions or adjust as needed
    # You can specify other extensions like '.txt', '.data', etc.
    if not (f.lower().endswith('.tsv') or f.lower().endswith('.txt')):
        print(f"Skipping non-tab-delimited file: {f}")
        continue

    print(f"\nProcessing file: {f}")
    try:
        # Key change: Specify sep='\t' for tab-delimited files
        file_df = pd.read_csv(input_file_path, sep='\t')
        print(file_df.columns)
        print(f"  Initial rows: {len(file_df)}")
    except Exception as e:
        print(f"❌ Error reading {f}: {e}. Skipping this file.")
        continue

    # Drop rows where 'Cigar' column is NaN.
    # This ensures that only rows with CIGAR information are processed.
    initial_rows_after_read = len(file_df)
    file_df.dropna(subset=['Cigar'], inplace=True)
    if len(file_df) < initial_rows_after_read:
        print(f"  Dropped {initial_rows_after_read - len(file_df)} rows due to NaN in 'Cigar' column.")
    print(f"  Rows after dropping NaN in 'Cigar': {len(file_df)}")

    if file_df.empty:
        print(f"  No valid CIGAR data in {f} after cleaning. Skipping.")
        continue

    # Apply the mutation extraction function.
    try:
        stat_df = mut_extraction(file_df)
        print(f"  Mutation extraction complete for {f}.")
    except Exception as e:
        print(f"❌ Error during mutation extraction for {f}: {e}. Skipping this file.")
        continue

    # Construct the output file path.
    # The original code named output files as "analyzed_" + input_filename.
    # You might want to keep the output as CSV, or change it to TSV as well.
    # For now, it remains CSV for output. If you want TSV output, change .to_csv to .to_csv(sep='\t')
    output_file_name = "analyzed_" + os.path.splitext(f)[0] + ".csv" # Keep output as .csv
    output_file_path = os.path.join(args.output_folder, output_file_name)

    # Save the processed DataFrame to a new CSV file.
    try:
        stat_df.to_csv(output_file_path, index=False) # index=False prevents writing DataFrame index.
        print(f"  Analyzed data saved to: {output_file_path}")
        processed_count += 1
    except Exception as e:
        print(f"❌ Error saving output for {f} to {output_file_path}: {e}")

print(f"\n✅ Script finished. Successfully processed {processed_count} files.")