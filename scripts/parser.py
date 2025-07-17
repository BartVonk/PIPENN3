from sys import argv

# to parse: parser.py [A] [B] [C]
# A: FASTA file name and location, e.g. data/my_fasta_file.fasta
# B: desired CSV file name and location, e.g. data/my_new_csv_file.csv
# C: (OPTIONAL) use -p to use prediction if present in FASTA
# Examples terminal command:
# python parser.py data/normal_fasta.fasta data/normal_csv.csv
# python parser.py data/prediction_fasta.fasta data/prediction_csv.csv -p
# python parser.py data/my_fasta_file.fasta data/my_new_csv_file.csv -p

def read_fasta(fastaFile, usePrediction):
    # function for parsing FASTA data
    # This script is for converting PiSITE(seq-insite) data, which actually uses PDB Chain-IDs
    headers, sequences, lengths, prediction = [["uniprot_id"], ["sequence"], ["Rlength"], ["p_interface"]]
    minRlength = 100
    maxRlength = 0

    with open(fastaFile, 'r') as file:
        headerLine = None
        sequenceLine = None
        predictionLine = None
        lengthLine = None
        for line in file:
            if not line:
                continue
            line = line.strip()
            if line.startswith('>'):
                # header
                headerLine = (line[1:])
            elif line.startswith('0') or line.startswith('1'):
                # prediction: read if on, otherwise zeroes of appropriate length
                predictionLine= (str('"' + ','.join(line) + '"' if usePrediction else '"' +','.join(len(line) * '0') + '"'))
            else:
                # sequence
                sequenceLine = (str('"' + ','.join(line) + '"'))
                lengthLine = (len(line))
                # keep track of min, max  of seq lengths
                if len(line) < minRlength:
                    minRlength = len(line)
                if len(line) > maxRlength:
                    maxRlength = len(line)
                
            #print(headerLine, sequenceLine, predictionLine, lengthLine)
            if all(var is not None for var in (headerLine, sequenceLine, predictionLine, lengthLine)):# and lengthLine >= 50:
                headers.append(headerLine)
                sequences.append(sequenceLine)
                prediction.append(predictionLine)
                lengths.append(lengthLine)
                
                headerLine = None
                sequenceLine = None
                predictionLine = None
                lengthLine = None
    
    # fill zeroes if prediction is empty
    if len(prediction) == 1:
        for i in range(len(sequences)-1):
            prediction.append(','.join(len(sequences[i+1])*'0'))

    return [headers, sequences, lengths, prediction, minRlength, maxRlength]

def calc_norm(minRlength, maxRlength, lengths):
    # this function calculates the normalized length for every sequence length in the list 'lengths'
    normalized = ["normalized_length"]
    for seqLength in lengths[1:]:
        # list starts one index later due to the fact the first index is the column header
        aaNormList = []
        for i in range(seqLength):
            aaNormList.append((seqLength - minRlength) / float(2050 - minRlength))
            #aaNormList.append((seqLength - minRlength) / float(maxRlength - minRlength))
        thing = ','.join([str(x) for x in aaNormList])
        use = '"' + thing + '"'
        normalized.append(use)
    return normalized

def write_csv(fastaData, csvFile):
    # function for rewriting FASTA data into csv format
    # all parsed information from the input file is found in 'fastaData', while 'csvFile' contains the name of the new csv file
    with open(csvFile, 'w') as file:
        for i in range(len(fastaData[0])):
            row = [str(fastaData[j][i]) for j in range(len(fastaData))]
            file.write(','.join(row) + '\n')

def main():
    # this line checks for the optional -p flag in the command
    if len(argv) < 3:
        print("Parser.py: Please be sure to submit an input fasta location and output csv location.")
        print("e.g.: python parser.py data/normal_fasta.fasta data/normal_csv.csv")
        exit()
    # sets variables based on command-line arguments
    fastaFile, csvFile = argv[1], argv[2]

    # usePrediction: True if -p is third arg, otherwise False
    usePrediction = len(argv) > 3 and str(argv[3]) == "-p"
    # collects all data from the input file
    print("collecting fasta data")
    fastaData = read_fasta(fastaFile, usePrediction)
    # collects the minimum and maximum sequence lengths
    minRlength, maxRlength = fastaData[4:]
    # clean the dataset
    fastaData = fastaData[:4]
    # adds column of normalized lengths
    fastaData.append(calc_norm(minRlength, maxRlength, fastaData[2]))
    # switches to order of the columns to comply with PiPENN input format
    # writes the collected and formatted data to csv
    print("writing fasta data")
    write_csv(fastaData, csvFile)

main()