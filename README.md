# Calculate SARS-CoV-2 lineage-defining consensus mutations on the basis of a covsonar db

This light script relies on:
* python3
* pandas
* argparse

Installation and running:

```bash
# create a covsonar database (https://github.com/rki-mf1/covsonar)
# export the database:
python3 path/to/covsonar/sonar.py match --db mydb > out.csv

# install the script:
git clone https://github.com/jonas-fuchs/covsonar_con_mut
cd covsonar_con_mut

# run the script:
python3 mutfreq.py out.csv > results.tabular

# optional arguments:
-l/--lineages
# define lineages to calculate the consensus mutations on.
# If none are given the whole db is used.
# Multiple lineages can be separated by a comma
-p/--profile
# calculate aa or dna consensus
-t/--threshold
# 0-1 frequency threshold to be considered
-c/--combine-lineages
# True or False, if True calculates one consensus profile for all given lineages.
# Ignored if no lineages are given
```
