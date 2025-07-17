import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the original CSV file
df = pd.read_csv(input_file)

# Select only the columns you need
columns_to_keep = ["uniprot_id", "sequence", "Rlength", "p_interface", "normalized_length"]
df_selected = df[columns_to_keep]

# Save to a new CSV file
df_selected.to_csv(output_file, index=False)
