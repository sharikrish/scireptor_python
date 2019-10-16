#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de  -adapted from Katharina Imkeller-

DESCRIPTION

From the IgBLAST output select the information on VDJ segments that will be written into the database.

The program goes through the IgBLAST output and looks for each query identifier,
which is identical with the seq_id in reads or sequences table. For each seq_id,
it writes the corresponding segments of 1. and 2. IgBLAST rank to the database.
The corresponding VDJ_ids from the library are also looked up and inserted.
Orientation and locus are updated in the updatetable, according to whether IgBLAST used the reverse complement
and what locus the 1. hit V segment has. Version: IGBLASTN 2.2.28+

1. Open IgBLAST output.
2. Set up the necessary DB queries: insert to VDJ table, select from VDJ_library and update in sequences/reads table.
3. Initialize variables needed to parse from IgBLAST output.
4. Go through each line of IgBLAST output and collect information.
"""

import argparse
from argparse import RawTextHelpFormatter
import bcelldb_init as binit

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument("-io", "--igblastoutput", required=True, help="IgBLAST output file.")
parser.add_argument("-t","--targettable", required=True,
                    help= "where VDJ info will be written. Can be reads_VDJ_segments or VDJ_\
                        segments according to where the sequences in the IgBLAST analysis came from.")
parser.add_argument("-ut","--updatetable", required=True, help="Table that will be updated according to what was found\
                        by IgBLAST (orientation, locus). Can be reads or sequences.")
args = parser.parse_args()

### 0. Logging information (TBD) # get configuration using bcelldb_init

conf = binit.get_config()

db = binit.connect()
cursor = db.cursor()

### 1. Open IgBLAST output file

try:
    in_igblast = open(args.igblastoutput,"r")
except FileNotFoundError:
    print ('\n\tBLAST output {} not found\n').format(args.igblastoutput)

# ### 2. Set up database handles for inserting to VDJ table and updating position

ins_VDJ_sth = ('INSERT IGNORE INTO {}.{} (seq_id, type, locus,'
                            ' igblast_rank, name, eval, score, VDJ_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s);').format(
    conf['database'],args.targettable)

# get segment_id and locus from the library
# the new igblast version does not output the locus anymore :(

sel_libr_sth = ('SELECT VDJ_id,'
                ' locus FROM {}.VDJ_library WHERE species_id=\"{}\" AND seg_name="%s";').format(
    conf["library"],conf["species"])

# update orientation and locus only
upd_posloc_sth = ('UPDATE {}.{} SET orient=%s, locus=%s, igblast_productive=%s WHERE seq_id=%s;').format(
    conf["database"],args.updatetable)

# count inserted seq_ids
count_ins = 0
# count total number of seq_id in blastoutput
count_total = 0

VDJ_id = ""
VDJ_locus = ""

print (conf["database"])
igblast_productive = "NULL"


for line in in_igblast:
    seq_orient = "F"
    count_total += 1

    # identify the query and extract id
    if "Query:" in line:
        query_id = line.strip().split(" ")[2]

    # identify start of rearrangement summary, parse following data
    if "rearrangement summary" in line:
        rea_fields = (next(in_igblast).strip().split("\t"))
        if rea_fields[-2] == "Yes":
            igblast_productive = 1
        elif rea_fields[-2] == "No":
            igblast_productive = 0

    # identify hit table and extract number of hits
    if "hits found" in line:
        count_V, count_D, count_J = 0, 0, 0
        n_hits = line.strip().split(" ")[1]

        # read "hit lines" until max number of hits reached
        for hit in range(0,int(n_hits)):
            fields=(next(in_igblast).split("\t"))
            VDJ_type = fields[0]
            seq_id = fields[1]
            VDJ_name = fields[2]
            evalue = fields[11]
            score = fields[12].strip()

            if ":" in VDJ_name: # needed for mouse GRCm38 database, where position information is part of the FASTA full name
                VDJ_name=VDJ_name.split(":")[0]

            if "reversed" in seq_id:
                seq_orient = "R"
            # collect 2 segments of each type
            if "V" in VDJ_type and count_V <= 1:
                count_V+=1
                # select id and locus from library
                cursor.execute(sel_libr_sth % VDJ_name)
                VDJ=cursor.fetchall()
                for id in VDJ:
                    VDJ_id=id[0]
                    VDJ_locus = id[1]
                    cursor.execute(ins_VDJ_sth, (query_id,VDJ_type,VDJ_locus,count_V,VDJ_name,evalue,score, VDJ_id))
                if count_V == 1:
                    # update orientation and locus only
                    cursor.execute(upd_posloc_sth,  (seq_orient,VDJ_locus,igblast_productive,query_id))

            elif "D" in VDJ_type and count_D <= 1:
                count_D += 1
                cursor.execute(sel_libr_sth % (VDJ_name))
                VDJ = cursor.fetchall()
                for id in VDJ:
                    VDJ_id = id[0]
                    VDJ_locus=id[1]
                    cursor.execute(ins_VDJ_sth, (query_id,VDJ_type,VDJ_locus,count_D,VDJ_name,evalue,score, VDJ_id))

            elif "J" in VDJ_type and count_J <= 1:
                count_J += 1
                cursor.execute(sel_libr_sth % (VDJ_name))
                VDJ = cursor.fetchall()
                for id in VDJ:
                    VDJ_id = id[0]
                    VDJ_locus = id[1]
                    cursor.execute(ins_VDJ_sth, (query_id,VDJ_type,VDJ_locus,count_J,VDJ_name,evalue,score, VDJ_id))

            elif VDJ_type not in ["V","D","J"]:
                print(
            "[todb_VDJ.py][WARNING] Encountered unknown segment type '{}' \
             in line: '{}'of input file. \n".format(count_total, VDJ_type))