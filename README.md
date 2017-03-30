
# Introduction
GRSR is a Tool for Deriving Genome Rearrangement Scenarios for Multiple Uni-chromosomal Genomes. This tool will do the following steps:
* Step 1. Run mugsy to get multiple sequence alignment results.
* Step 2. Get a permutation for each strain with some parameters.
* Step 3. Generate pairwise genome rearrangement scenarios and find repeats at the breakpoints of each rearrangement events.

The package works under Linux system.
* Folder scr contains code developed by me (Author: Dan WANG, danwang5-c@my.cityu.edu.hk).
* Folder Musgy contains code by Angiuoli SV and Salzberg SL. The code was downloaded from http://mugsy.sourceforge.net/
* Folder GRIMM-synteny and GRIMM contains code by Glenn Tesler. The code was downloaded from http://grimm.ucsd.edu/DIST/
* Folder blast-2.2.26 contains code by NCBI. The code was downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.26/ 

# Example  
Before running the code, please put all the five folders (scr,Musgy,blast-2.2.26,GRIMM-synteny and GRIMM) in the same directory.
Here is an example that you can run the three steps separately. All the inputs and outputs of each step are in the [~/src/Example](src/Example) directory. Step 1 may take a long time. For example, for aligning 25 bacterial chromosomes, it will take almost 1 day. But for aligning the 3 bacterial chromosomes in our example, it will only take several minutes. Step 2 and Step 3 are fast and normally will only cost a few minutes.
## Step 1. Run mugsy to get multiple sequence alignment results:
Commands for Step 1:
```
cd ~/Mugsy/mugsy_x86-64-v1r2.3
mugsy --directory ~/scr/Example/1.MSA_results --prefix MSA_Result ~/scr/Example/Genomes/*.fna
```
where
* [~/scr/Example/1.MSA_results](scr/Example/1.MSA_results): the direcotry for output file
* [MSA_Result](scr/Example/1.MSA_results/MSA_Result.maf): the prefix of the output file is MSA_Result
* [~/scr/Example/Genomes/*.fna](scr/Example/Genomes): all the input chromosomes which are in the fna format.(For the other input formats of mugsy, you can refer to the mugsy's website.)

This step outputs
 * [~/scr/Example/1.MSA_results/MSA_Result.maf](scr/Example/1.MSA_results/MSA_Result.maf): multiple sequence alignment results

## Step 2. Get a permutation for each strain with some parameters.
Commands for Step 2:
```
cd ~/scr
sh getScaffolds.sh 500 3000
```
where 
* 500: is the minimum block size
* 3000: is the maximum gap threshold

This step outputs
* [~/scr/Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt](scr/Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt): the scaffolds for each strain. (Each integer in the scaffolds stands for a core-genome block whose length is larger than 500. And consecutive blocks in all the strains are merged into one block.)
* [~/scr/Example/2.Scaffolds/anchors_500_3000_c/blocks.txt](scr/Example/2.Scaffolds/anchors_500_3000_c/blocks.txt): keeps the position of each block in the mgr_macro.txt file
* [~/scr/Example/2.Scaffolds/cordis.maf](scr/Example/2.Scaffolds/cordis.maf): the multiple sequence alignment result without the strain-specific segments
* [~/scr/Example/2.Scaffolds/core_coords.txt](scr/Example/2.Scaffolds/core_coords.txt): is the positions of core-genome blocks in each strain without filtering short blocks and merging consecutive blocks.
 
## Step 3. Generate pairwise genome rearrangement scenarios and find repeats at the breakpoints of each rearrangement events.
- Firstly, merget the consecutive blocks in both the two strains;
- Secondly, detect independent block-interchange and transposition events and find repeats at the breakpoints, the result is written to the reportA.txt;
- Thirdly, calculate the reversals between the two genomes by using grimm after removing the independent block-interchanges and transpositions in step 2;
- Lastly, check whether a pair of inverted repeats exist at the two ends of an reversal and store the results in reportA.txt

Commands for Step 3:
```
cd ~/scr
sh reptA.sh ./Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt ./Example/2.Scaffolds/anchors_500_3000_c/blocks.txt .Example/2.Scaffolds/cordis.maf ./Example/Genomes ./Example/3.Repreport 100 90 5000
```
where
* [./Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt](scr/Example/2.Scaffolds/anchors_500_3000_c/mgr_macro.txt): the scaffolds for each strain.
* [./Example/2.Scaffolds/anchors_500_3000_c/blocks.txt](scr/Example/2.Scaffolds/anchors_500_3000_c/blocks.txt): the position of each block on the strains' genomes
* [./Example/2.Scaffolds/cordis.maf](scr/Example/2.Scaffolds/cordis.maf): the multiple sequence alignment result without the strain-specific segments
* [./Example/Genomes](scr/Example/Genomes): the directory where input chromosome sequences are stored.
* [./Example/3.Repreport](scr/Example/3.Repreport): the output files directory 
* 100: The minimum length of repeats. if the length of a repeat is less than 100, we will igore it.
* 90: The minimum similarity between repeats. If the similarity between two repeats are smaller than 90, we will ignore it.
* 5000: The range we search for repeats at the breakpoints of reversals.

This step outputs
[~/scr/Example/3.Repreport/reportA.txt](scr/Example/3.Repreport/reportA.txt): provides the rearrangement scenario between any pair of strains and whether repeats exists at the two ends of an reversal. For example: The following paragraph means that to transform Strain 2 (source) to Strain 3 (destination), there are a total of 11 rearrangement steps (1 transposition, 3 inverted transpositions and 7 reversals) from Strain 2 to 3. The transposition is Block 20 and no repeats are found at its breakpoints in both Strain 2 and 3. The three inverted transposition is Block 2, 5 and 24 and no repeats are found at their breakpoints. Aftering removing the blocks involved in transpositions and block interchange, there are 7 reversals from Starin 2 to 3. The first reversal is to inverse Block 23 trough Block 25 in strain 2 and no pair of inverted repeats (A/-A) are found at the two ends of this reversal in both source and destination's genomes. The third reversal is to inverse Block 7 trough Block 21 in strain 2. At the two ends of this reversal, a pair of IRs are found in Strain 2 (source), the length of this IR(A/-A) is 5001 bp and the similarity betWeen +A and -A is 99%. Also, a pair of IRs are found in Strain 3 (destination), the length of this IR(A/-A) is 5000 bp and the similarity betWeen +A and -A is 99%.
-------------------------------------------------------------------------------------
	>From Genome 2 to 3, total rearrangment step(s): 11
	*Transposition: 1 step(s):
	Step 1: 20 through 20 Transposition:Source:N Destination:N
	*Inverted Transposition: 3 step(s):
	Step 1: 2 through 2 Inverted Transposition: Source:N Destination:N
	Step 2: 5 through 5 Inverted Transposition: Source:N Destination:N
	Step 3: 24 through 24 Inverted Transposition: Source:N Destination:N
	*Reversal: 7 step(s):
	Step 1: 23 through 25 Reversal: Source:N Destination:N
	Step 2: -8 through -23 Reversal: Source:N Destination:N
	Step 3: 7 through 21 Reversal: Source:Found,Length = 5001,Simlarity = 99% Destination:Found,Length = 5000,Simlarity = 99%
	Step 4: -16 through 22 Reversal: Source:N Destination:N
	Step 5: 9 through 25 Reversal: Source:N Destination:N
	Step 6: -28 through 8 Reversal: Source:N Destination:N
	Step 7: -22 through 18 Reversal: Source:N Destination:N
------------------------------------------------------------------------------------
