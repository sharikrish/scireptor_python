#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de  -adapted from Katharina Imkeller-

DESCRIPTION
Open a BLAST output file. Parse the necessary information and write to database.
The identifiers in the BLAST output need to be the same as the seq_id in reads or sequences table.

1. Prepare DB statements: insert into constant table, select id from constant library.
2. Open BLAST output.
3. Parse information from the output.
"""

import argparse
from argparse import RawTextHelpFormatter
import bcelldb_init as binit

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument("-bo", "--blastoutput", help="blast output file")
parser.add_argument("-t","--targettable", help="total count of sequences")
args = parser.parse_args()

### 0. Logging information (TBD) # get configuration using bcelldb_init

conf = binit.get_config()

db = binit.connect()
cursor = db.cursor()

### 1. Prepare statement to get constant_id from library and insert constant information

# statement to find out which is the corresponding constant_id in the library
# name in the constant_library needs to be unique!!!


### 2. Open BLAST output

try:
    blast = open(args.blastoutput,"r")
except FileNotFoundError:
    print (f'\nBLAST output {args.blastoutput} not found')


### 3. Parse information

# count inserted seq_ids
count_ins = 0
# count total number of seq_id in blastoutput
count_total = 0

for line in blast:
    # increase total counting
    count_total +=1

    (seq_id, hit_name_big, percid, length,
     mismatches, gapopens, readstart, readend,
     conststart, constend, evalue, score) = line.strip().split("\t")

    # find out which constant segement was found

    # in the mouse database, there is some extra information in the identifier
    if ':' in hit_name_big:
        hit_name = (hit_name_big.split(':')[0])

    else:
        hit_name = hit_name_big

    # statement to find out which is the corresponding constant_id in the library
    # name in the constant_library needs to be unique!!!

    # select the corresponding constant_id from the database
    sel_statem = f'SELECT constant_id FROM {conf["library"]}.constant_library ' \
                 f' WHERE species_id="{conf["species"]}" AND name="{hit_name}"'
    cursor.execute(sel_statem)
    constant_id = cursor.fetchall()
    for id in constant_id:
        constant_id=id[0]
    # insert into database
        ins_statem = f'INSERT IGNORE INTO {conf["database"]}.{args.targettable} (seq_id,' \
        f' name, percid, length, gapopens, readstart, readend, eval,' \
        f' score, constant_id) VALUES ("{seq_id}","{hit_name}","{percid}","{length}",' \
        f'"{gapopens}","{readstart}","{readend}","{evalue}","{score}","{constant_id}")'
        ins_bool=cursor.execute(ins_statem)

    # update inserted count
        count_ins +=ins_bool


print (f"\n---------\ntotal {count_total} sequences processed\n")
print (f"{count_ins} inserted\n---------\n\n")

blast.close()