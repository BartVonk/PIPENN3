import os
from sys import argv

def fetchFiles(location):
    # This function generates a list of all files in the data/pisite/all folder
    fileList = []
    for root, _, files in os.walk(location):
        for file in files:
            fileList.append(os.path.join(root, file))
            #return(fileList)

    return(fileList)

def scanPiSITE(file):
    # This function parses a single pisite file for the PDBID, CHAINID, Sequence, and presence of binding partners

    # initiate / reset values
    sequence = []#3rd column
    any_pred = []#4th column
    PDBID = None
    CHAIN = None
    startParsing = False

    # opens the provided pisite filename
    with open(file, "r") as f:
        for line in f:
            # reads the pisite file and scans for specific information markers
            if line.startswith('#PDBID:'):
                # PDBID marker
                PDBID = line[8:].strip()
            elif line.startswith('#CHAIN:'):
                # CHAIN marker
                CHAIN = line[8:].strip()
                PDBID = PDBID+CHAIN
            elif line.startswith("#residue_no insertion_code residue_type number_of_binding_partners binding_partners"):
                # marker for start of sequence and binding partner columns
                startParsing = True
                continue

            if startParsing:
                parts = line.split()
                if len(parts) >= 4:
                    # extracts sequence (second column when parsing has started)
                    sequence.append(parts[2])
                    # extracts binding partner string and transforms it to binary, 0 if absent, 1 if one or more is present
                    any_pred.append('0' if parts[3] == '0' else '1')

    return(PDBID.rstrip(), "".join(sequence), "".join(any_pred))

def main():
    pisiteLocation = 0

    if len(argv) < 3:
        print("python Parser.py pisite/all/ data/out.fasta")
        exit()

    # collects pisite filenames
    fileList = fetchFiles(argv[1])
    outLoc = argv[2]
    total_files = len(fileList)

    # Define thresholds at every 10%
    progress_thresholds = {i for i in range(10, 101, 10)}
    last_reported = 0

    with open(outLoc, "w") as outputFile:
        for idx, file in enumerate(fileList):
            PDBID, sequence, any_pred = scanPiSITE(file)
            outputFile.write(f">{PDBID}\n{sequence}\n{any_pred}\n")

            percent_done = int(((idx + 1) / total_files) * 100)
            if percent_done >= last_reported + 10 and percent_done in progress_thresholds:
                print(f"Progress: {percent_done}% complete")
                last_reported = percent_done



main()

