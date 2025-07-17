from sys import argv

# This script is made to remove all data entries from a data that match a list of IDs (e.g remove all IDs in zk448 from our dataset)
# python filter_protein_entries.py ids.txt entries.fasta output.fasta

def read_ids(file):
    """Read IDs from a file and return a set of IDs."""
    with open(file, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def filter_entries(ids_file, entries_file, output_file):
    """Filter entries based on IDs and save the remaining ones."""
    idsToRemove = read_ids(ids_file)
    #print(idsToRemove)
    
    with open(entries_file, 'r') as f, open(output_file, 'w') as out:
        keepEntry = False
        entryLines = []
        
        for line in f:
            if line.startswith('>'):
                if entryLines and keepEntry:
                    out.writelines(entryLines)  # Write previous entry if it should be kept
                
                entry_id = line[1:].strip()  # Remove '>' and strip whitespace
                keepEntry = entry_id not in idsToRemove
                entryLines = [line] if keepEntry else []
            elif keepEntry:
                entryLines.append(line)
        
        if entryLines and keepEntry:  # Write the last entry if it should be kept
            out.writelines(entryLines)

if __name__ == "__main__":
    if len(argv) != 4:
        print("Usage: python script.py <ids_file> <entries_file> <output_file>")
        exit()
    
    ids_file, entries_file, output_file = argv[1], argv[2], argv[3]
    filter_entries(ids_file, entries_file, output_file)