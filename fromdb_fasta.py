#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de  -adapted from Katharina Imkeller-

DESCRIPTION

Get a fasta file of sequences in the sourcetable that do not have any correspondence in the targettable. Modulo and rest are used to distribute all the sequences to smaller files that can then be analysed by different processes.
Needed as a step before executing further analysis.
The fasta file will be used as source for IgBLAST, BLAST, RazerS and MUSCLE. The fasta file is written into the output subrepository.

1.	Parse segmentation variables from the outfile variable.
2.	Select all entries of sourcetable that do not have corresponding entry in targettable.
3.	Write seq_id and seq of selected rows into fasta file.

- select statement for sequences
- total number of sequences that were selected

"""

import argparse
from argparse import RawTextHelpFormatter
import re
import bcelldb_init as binit


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument("-s", "--sourcetable", help="Source table")
parser.add_argument("-t", "--targettable", help="Target table")
parser.add_argument("-f", "--outfile", help="Output file")
args = parser.parse_args()

source_table = args.sourcetable
target_table = args.targettable
outfile = args.outfile

### 0. Logging information (TBD) # get configuration using bcelldb_init

conf = binit.get_config()

database_name = conf['database']

db = binit.connect()
cursor = db.cursor()

### 1. Parse segmentation variables from the outfile variable.

regex_pattern = "\d+_\d+" #works with this regex, check differences
regex = re.compile(regex_pattern)
mod, rest = (regex.findall(outfile)[0].split("_"))

### 2. Get all entries from sourcetable that do not have correspondence in targettable

# prepare statement to get sequences

statement=f"SELECT input.seq_id, input.seq \
  FROM {database_name}.{source_table} AS input \
  LEFT JOIN {database_name}.{target_table} AS output ON output.seq_id=input.seq_id \
  WHERE output.seq_id IS NULL"

if int(mod) > 1:
    statement = statement + " and MOD(input.seq_id,{}) = {} ".format(mod,rest)

# log the select statement
print(f"\n\nselect statement: {statement}\n\n")

# get sequences

seq_count = 0
cursor.execute(statement)

### 3. Write seq_id and seq of selected rows into fasta file.

with open (outfile, "w") as fasta:
    for line in cursor.fetchall():
        seq_count+=1
        # print(">{}\n{}\n\n".format(line[0],line[1]))
        fasta.write(">{}\n{}\n\n".format(line[0],line[1]))

print(seq_count)
