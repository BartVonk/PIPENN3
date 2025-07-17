import sys
import csv
import os
import random

def split_csv(input_file, output_file_70, output_file_30, split_ratio=0.95):
    with open(input_file, 'r', newline='', encoding='utf-8') as infile:
        reader = list(csv.reader(infile))  # Read all rows into a list
        
        if not reader:
            print("Error: The CSV file is empty.")
            return
        
        header, rows = reader[0], reader[1:]
        random.shuffle(rows)  # Shuffle rows to ensure randomness
        
        split_index = int(len(rows) * split_ratio)
        rows_70, rows_30 = rows[:split_index], rows[split_index:]
        
        # Write the 70% split file
        with open(output_file_70, 'w', newline='', encoding='utf-8') as out70:
            writer = csv.writer(out70)
            writer.writerow(header)
            writer.writerows(rows_70)
        
        # Write the 30% split file
        with open(output_file_30, 'w', newline='', encoding='utf-8') as out30:
            writer = csv.writer(out30)
            writer.writerow(header)
            writer.writerows(rows_30)
        
    print(f"Successfully split '{input_file}' into:")
    print(f"  - {output_file_70} ({len(rows_70)} rows)")
    print(f"  - {output_file_30} ({len(rows_30)} rows)")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python myscript.py myfile.csv")
        sys.exit(1)
    
    input_filename = sys.argv[1]
    if not os.path.exists(input_filename):
        print(f"Error: File '{input_filename}' not found.")
        sys.exit(1)
    
    base_name = os.path.splitext(input_filename)[0]  # Extract filename without extension
    output_training = f"{base_name}_training.csv"
    output_testing = f"{base_name}_testing.csv"
    
    split_csv(input_filename, output_training, output_testing)
