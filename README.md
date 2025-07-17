# PIPENN3
this pipeline is designed to convert an unzipped release from pisite, whose file structure is located inside a /all/ folder, into an prepared csv file to be used by PIPENN.
to achieve this the pipeline requires two outside elements:
--> the path to the pisite /all/ folder.
--> a single textfile containing PDB ids that are not allowed in the final dataset, the pipeline takes homologs into account.

to start pipeline run the following code from the desired location: . ../../pipeline/GenerateData.sh ../../data/pisite/all/ BioDLIDs_PDB.txt

inside the pipeline there are two variables that can be changes from the default set value:
--> the minimum sequence length
--> the minimum sequence identity

the minimum sequence length is defined as the number at the end of the line calling 'generalFilter.py' (line 27 at the time of writing), this variable is 50 by default; sequences smaller than 50 AA's will be deleted, there is no maximum lenght variable.

the minimum sequence identity is used at two separate times, initially for the cluster development and later on when fetching the representative sequences. this value is currently set to 0.25 (lines 37 & 47)
