#!/usr/bin/python

####################################
# File name: todb_igblast_align.py#
# Author: Srilakshmy               #
# Created Date: August 26th 2019   #
####################################

import argparse
import bcelldb_init as binit
import re
##usage : todb_igblast_align.py [-h] -io <igblastoutput> -dir <directory_for_output>

##get configuartion for scireptor

config_scireptor=binit.get_config()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-io', '--igblastoutput', required=True, help="define input file")
    parser.add_argument('-dir', '--directory', required=True, help="define output file")
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    igblastopt = args.igblastoutput
    outdir = args.directory
###Prepare database for insertion of sequence
### Logging and database init
    conn = binit.connect()
##create a cursor object
    cursor = conn.cursor()
###prepare insert command
    ins_aln_statement = 'INSERT IGNORE INTO {}.igblast_alignment (seq_id, query_start, germline_start, query_seq, germline_seq) VALUES (%s,%s,%s,%s,%s);'.format(config_scireptor['database'])
    # go through file
    with open(igblastopt) as f:
        lines=f.readlines()
        for i, line in enumerate(lines):
            if "# Hit table" in line:
                fields=lines[i+1].split(",")
                if not (re.match("query id", fields[0]) and
                    re.match("q. start", fields[6]) and
                    re.match("s. start",fields[8]) and
                    re.match("query seq",fields[12]) and
                re.match("subject seq", fields[13])):
                    print("Fields in the hit table of IgBLAST output changed? Cannot parse.")

                first_hit=lines[i+3].split("\t")

                # Read data elements. ATTENTION: Output elements are shifted to the right by 1, due to
            # an additional letter in the raw output, which indicates the segment type.
                query_id = first_hit[1]
                query_start = first_hit[7]
                germline_start = first_hit[9]
                query_seq = first_hit[13]
                germline_seq = first_hit[14]
                val = (query_id, query_start, germline_start, query_seq, germline_seq)
                # insert into database
                ####REMEMBER need todo santization###
                cursor.execute(ins_aln_statement,val)

