# Find ORFs within a virus genome
It allows you to find orfs in fasta sequences of viral genomes.

If a reference is provided as genebank file, the script uses the min length of the orfs within the genebank file and determines automatically if the orfs should be overlapping or not.

## This script relies on:
* python3
* pandas
* SeqIO


## Which orfs can the script find:

```bash
Types:  
-----[M-------------*]-------- complete  
------M--[M---------*]-------- complete_internal  
-----[M----------------------- 5_partial  
------M--[M------------------- 5_partial_internal  
[--------------------*]------- 3_partial  
[----------------------------] 5_3_partial  
------*]-------------[M------- circular  
```

## How the no overlap algorithm works:

```bash
no overlap algorithm:  
frame 1: -[M------*]-------[M--*]---------[M-----  
frame 2: -------[M------*]---------[M---*]--------  
frame 3: [M---*]-----[M----------*]----------[M---  

results: [M---*][M------*]--[M--*]-[M---*]-[M-----  
frame:    3      2           1      2       1  
```

## Installation and running:  

```bash

# install the script:
git clone https://github.com/jonas-fuchs/viral_orf_finder
cd viral_orf_finder

# run the script:
python3 orf_finder.py infile.fasta

# optional arguments:
-r/--reference
# reference genebank file
-m/--min-length
# min length of the orfs to find
-i/--internal
# True/False if script should search for internal orfs
-c/--circular
# True/False if script should search for circular orfs
-p/--partial
# True/False if script should search for partial orfs
# with no START and/or STOP
-n/--no-overlap
# True/False if script should consider only orfs that do not overlap
-s/--strands
# + and/or - (+ = positive strand, - = negative strand)
```

### Side note:

The script has set the following codons:
start_codons = ["ATG"]
stop_codons = ["TAG", "TGA", "TAA"]

If you want to edit this you can do this directly by setting codons in the find_orfs function within the script.
