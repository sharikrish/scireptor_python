#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de  -adapted from Katharina Imkeller-

DESCRIPTION 

Extract all the tags found by RazerS and write them to the database. Calculate number of insertions, deletions and mutations in the tag.

Steps:
1. Prepare database to insert tag information
2. Open RazerS output
3. Parse tag information and write to database

LOGGED INFORMATION

- count of inserted tags
"""

import argparse
from argparse import RawTextHelpFormatter
import bcelldb_init as binit


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument("-ro", "--razerout", help="Razerout file")
args = parser.parse_args()
razersout = args.razerout


# get database handle

db = binit.connect()
cursor = db.cursor()

conf = binit.get_config()

### 1. Prepare database to insert tag information


### 2. Open RazerS output

try:
    razer = open(razersout, "r")
except IOError:
    print(f"opening razers output {razersout} failed")

count =0

for line, hash_read, genome in zip(razer,razer,razer):

    count+=1

    # print (line,read, genome)

    tag,  tag_start, tag_end, orient, read, read_start, read_end, percid = (line.strip().split("\t"))
    tag_seq = hash_read.split(":")[1].strip()
    read_seq = genome.split(":")[1].strip()
    if percid != 100:
        nins = tag_seq.count("-")
        # print (nins)
        ndel = read_seq.count("-")
        nmut=len(tag_seq)*(1-float(percid)/100)-nins -ndel

# insert the tag information
    insert_tags_statement = f'INSERT IGNORE INTO {conf["database"]}.reads_tags (' \
        f'seq_id, percid, direction, insertion, deletion, replacement, start, end,' \
        f' tag_id)VALUES ("{read}","{percid}","{orient}","{nins}","{ndel}","{nmut}","{read_start}","{read_end}","{tag}")'

    cursor.execute(insert_tags_statement)

razer.close()

print(f"\n\nA total number of {count} tags were inserted.\n\n")