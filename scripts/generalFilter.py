import sys

def filter_fasta(input_file, output_file, min_length, max_length):
    # Open the input and output files
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Initialize variables
        entry_id = None
        sequence = None
        binary_string = None
        
        # Iterate through the lines in the input file
        for line in infile:
            line = line.strip()
            
            # If the line starts with '>', it means we are reading an ID line
            if line.startswith(">"):
                if entry_id and sequence and binary_string:
                    # Check if the current entry meets the filtering criteria
                    if len(sequence) >= min_length and len(sequence) <= max_length and binary_string != '0' * len(binary_string) and binary_string != '1' * len(binary_string):
                        # If so, write it to the output file
                        interface_perc = (binary_string.count('1') / len(binary_string)) * 100 if binary_string else 0
                        if len(sequence) >50 and interface_perc >90:
                            pass
                        else:
                            if 'U' not in sequence:
                                outfile.write(f"{entry_id}\n{sequence}\n{binary_string}\n")
                            else:
                                pass
                
                # Start a new entry
                entry_id = line
                sequence = ""
                binary_string = ""
            
            elif line.isalpha():
                # Add to the sequence
                sequence += line
            elif line.isdigit():
                # Add to the binary string (should be a string of 0s and 1s)
                binary_string += line
        
        # After the loop, we need to check the last entry
        if entry_id and sequence and binary_string:
            if len(sequence) >= min_length and binary_string != '0' * len(binary_string):
                outfile.write(f"{entry_id}\n{sequence}\n{binary_string}\n")

if __name__ == "__main__":
    # Ensure the script is called with the right number of arguments
    if len(sys.argv) != 5:
        print("Usage: python generalFilter.py <input.fasta> <output.fasta> <min_length> <max_length>")
        sys.exit(1)

    # Parse the arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_length = int(sys.argv[3])
    max_length = int(sys.argv[4])

    # Call the filter function
    filter_fasta(input_file, output_file, min_length, max_length)
