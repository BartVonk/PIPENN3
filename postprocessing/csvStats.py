from sys import argv
from re import findall
import csv
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# python csvStats.py ../PiSITE/Seq-InSITE.csv

def count_consecutive_ones(INTERFACE):
    """Counts occurrences of consecutive '1's and finds the longest chain of '1's."""
    INTERFACE = INTERFACE.replace(",", "")
    ones_groups = findall(r'1111+', INTERFACE)  # Find all groups of consecutive '1's, at least size 4
    num_consecutive_ones = len(ones_groups)  # Number of groups of '1's
    max_chain = max(map(len, findall(r'1+', INTERFACE)), default=0)  # Longest stretch of '1's
    return num_consecutive_ones, max_chain

def getStats(input_file):
    statsDict = {"minLength": 0,            # Min Length
                 "maxLength": 0,            # Max Length
                 "meanLength": 0,           # Sum Length / nr of Seqs
                 "yesInterface": 0,         # Interface Ratio
                 "noInterface": 0,          # Interface Ratio
                 "avgMaxConsInterface": 0,  # average (file) maximum consecutive interface (seq)
                 "nrConsInterface": 0,      # Total
                 "InterFacePerAA": 0}       # 
    
    amino_acids_interface = {
    'A': 0,  # Alanine
    'R': 0,  # Arginine
    'N': 0,  # Asparagine
    'D': 0,  # Aspartic acid
    'C': 0,  # Cysteine
    'E': 0,  # Glutamic acid
    'Q': 0,  # Glutamine
    'G': 0,  # Glycine
    'H': 0,  # Histidine
    'I': 0,  # Isoleucine
    'L': 0,  # Leucine
    'K': 0,  # Lysine
    'M': 0,  # Methionine
    'F': 0,  # Phenylalanine
    'P': 0,  # Proline
    'S': 0,  # Serine
    'T': 0,  # Threonine
    'W': 0,  # Tryptophan
    'Y': 0,  # Tyrosine
    'V': 0,   # Valine
    'X': 0
}
    amino_acids_general = {
    'A': 0,  # Alanine
    'R': 0,  # Arginine
    'N': 0,  # Asparagine
    'D': 0,  # Aspartic acid
    'C': 0,  # Cysteine
    'E': 0,  # Glutamic acid
    'Q': 0,  # Glutamine
    'G': 0,  # Glycine
    'H': 0,  # Histidine
    'I': 0,  # Isoleucine
    'L': 0,  # Leucine
    'K': 0,  # Lysine
    'M': 0,  # Methionine
    'F': 0,  # Phenylalanine
    'P': 0,  # Proline
    'S': 0,  # Serine
    'T': 0,  # Threonine
    'W': 0,  # Tryptophan
    'Y': 0,  # Tyrosine
    'V': 0,   # Valine
    'X': 0
}

    listInterfacePercentage = []
    listLengths = []
    with open(input_file, 'r') as file:
        next(file) # skips header
        minLength = 10000
        maxLength = 0
        sumLength = 0
        nrOfSeqs = 0
        yesInterface = 0
        noInterface = 0
        
        avgMaxConsInterface = 0
        minMaxConsInterface = 10000
        maxMaxConsInterface = 0
        nrConsInterface = 0
        poolLargest = 0
        for line in file:
            parsed_line = next(csv.reader([line]))

            ID = parsed_line[0]  # First value
            SEQ = parsed_line[1]  # Large sequence of letters
            LENGTH = parsed_line[2]  # Integer value
            INTERFACE = parsed_line[3]  # Large sequence of binary digits
            NORM_LENGTH = parsed_line[4]  # Floating point value

            # minLength, maxLength
            if int(LENGTH) < minLength:
                minLength = int(LENGTH)
            if int(LENGTH) > maxLength:
                maxLength = int(LENGTH)
            sumLength += int(LENGTH)
            nrOfSeqs += 1

            # interface
            noInterface += INTERFACE.count('0')
            yesInterface += INTERFACE.count('1')

            percentage = (round(INTERFACE.count('1')/int(LENGTH), 5)*100)
            if int(percentage) <= 1 :
                print(ID, percentage, LENGTH)

            listInterfacePercentage.append(percentage)
            listLengths.append(int(LENGTH))

            # consecutive interface
            num_consecutive_ones, max_chain_ones = count_consecutive_ones(INTERFACE)
            poolLargest += max_chain_ones
            nrConsInterface += num_consecutive_ones
            if int(max_chain_ones) < minMaxConsInterface:
                minMaxConsInterface = int(max_chain_ones)
            if int(max_chain_ones) > maxMaxConsInterface:
                maxMaxConsInterface = int(max_chain_ones)

            filtered_chars = [char for char, bit in zip(SEQ, INTERFACE) if bit == '1']
            #print(filtered_chars)
            for char in filtered_chars:
                if char in amino_acids_interface:  # Ensure the character is a valid amino acid
                    amino_acids_interface[char] += 1
                elif char not in amino_acids_interface.keys():
                    print("Amino Acid not found:", char)

            for char in SEQ:
                if char in amino_acids_general.keys():
                    amino_acids_general[char] += 1

        meanLength = sumLength / nrOfSeqs
        avgMaxConsInterface = poolLargest / nrOfSeqs

        # STATS 
        print("#csvStats of ", input_file, "\n")
        print("#Diff yInt & dict:\t", yesInterface-sum(amino_acids_interface.values()))
        print("Nr Seqs:\t\t", nrOfSeqs)                                                 # Total number of amino acids
        print("MeanLen:\t\t", round(meanLength, 3))                                     # average sequence length
        print("MinLen:\t\t\t", round(minLength, 3)) 
        print("MaxLen:\t\t\t", round(maxLength, 3)) 
        print("Sum length AA:\t\t", sumLength)                                            # tot sequence length
        print("yInt + nInt\t\t", (yesInterface + noInterface))                          # sum yInt + nInt (should be equal to sumLength)
        print("yInt:\t\t\t", yesInterface)                                              # L-> number of occurances of 1 in interface string
        print("nInt:\t\t\t", noInterface)                                               # L--> number of occurances of 0 in interface string
        print("Ratio Int/total:\t", str(round(yesInterface/sumLength, 5)*100)+ "%")         # Ratio interface / total
        print("Log-Odds Int/noInt:\t", str(round(math.log(noInterface / yesInterface), 5)))   # log-odds interface / no interface
        print("Nr ConsInt:\t\t", nrConsInterface)                                       # number of consecutive interfaces
        print("AvgMaxConsInt:\t\t", round(avgMaxConsInterface, 3))                      # average largest consecutive interface size
        print("MinMaxConsInt:\t\t", round(minMaxConsInterface, 3))
        print("MaxMaxConsInt:\t\t", round(maxMaxConsInterface, 3))
        print("Nr AA in int: sum dict\t", sum(amino_acids_interface.values()))                           # nr of interfaces, calculated from AA dict
         
        #print('2nd Int/total Ratio:\t', str(round(sum(amino_acids_interface.values())/sumLength, 9)*100)+"%")  # 2nd check for ratio since AA dict calc is off
        print('')                                                                       # 
        print("Breakdown of AA counts in interface:")#, amino_acids_interface)                          # breakdown of amino acid distribution in interaces
        for amino, count in amino_acids_interface.items():
            print(f"{amino}: {count}")

        print('')                                                                       # 
        print("Breakdown of AA counts in general:")#, amino_acids_interface)                          # breakdown of amino acid distribution in interaces
        for amino, count in amino_acids_general.items():
            print(f"{amino}: {count}")
        
        #array_data = np.array(list(listInterfacePercentage))  # forces materialization
        #array_data = np.array(list(listLengths))  # forces materialization
        #sns.kdeplot(array_data, bw_method=0.05)
        #sns.set_style('whitegrid')
        #sns.kdeplot(np.array(listInterfacePercentage), bw_method=0.05)
        #sns.kdeplot(np.array(listLengths), bw_method=0.05)

        #plt.savefig("intfacePercentageDensityPlot.png")  # Save the plot as a PNG file
        #plt.savefig("lengthDensityPlot.png")  # Save the plot as a PNG file
        #plt.close()  # Optional: close the plot window if running in a loop or script
        #listInterfacePercentage.sort()
        #print(listInterfacePercentage)

        violin_df = pd.DataFrame({
            'Length': listLengths,
            'Interface %': listInterfacePercentage
        })
        violin_df['Length_bin'] = pd.cut(violin_df['Length'], bins=10)
        plt.figure(figsize=(12, 6))
        violin_df['bin_label'] = violin_df['Length_bin'].apply(lambda x: f"{int(x.left)}–{int(x.right)}")
        sns.violinplot(x='bin_label', y='Interface %', data=violin_df, inner='quartile')
        plt.xticks(rotation=45)
        plt.xlabel('Length Bin')
        plt.ylabel('Interface %')
        plt.title('Interface Percentage Distribution Across Length Bins')
        plt.tight_layout()
        plt.savefig("violin_bins_interface_length.png")
        plt.close()



if __name__ == "__main__":
    if len(argv) != 2:
        print("Usage: python script.py <input_csv>")
        exit()
    
    input_file = argv[1]
    getStats(input_file)