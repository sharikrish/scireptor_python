#!/usr/bin/python3

####################################
# File name: fromdb_consensus_fasta.py#
# Author: Srilakshmy and Francisco  #
# Created Date: August 27th 2019   #
####################################

import sys, argparse, subprocess, os
from Bio.Seq import Seq
import bcelldb_init as binit
##usage -> fromdb_consensus_fasta.pl -p <output_path> [-h <help>]

##get configuartion for scireptor
conf = binit.get_config()
db = binit.connect()
n_consensus=conf['n_consensus']
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--output_path', required=True, help="define input file")
    args = parser.parse_args()

    output = args.output_path

### 1. Prepare database for insertion of sequence
    #conn = dbconnect.connect()
##create a cursor object
    cursor = db.cursor()
### Functions
    def revcomp(dna):
        seq_dna=Seq(dna)
        return seq_dna.reverse_complement()

### 1. Prepare select statements
# get all the consensus_ids, that do not have a matching sequence
    sel_consensus_sth = (f"SELECT consensus_id \
  FROM {conf['database']}.consensus_stats \
  WHERE sequences_seq_id IS NULL")

    cursor.execute(sel_consensus_sth)

### 2. Loop over all consensi
# counters for consensi, depending on wherther they have enough reads or not
    total_count = 0           # total number of selected consensi
    cons_count =0             #consensi that have more then threshold nb of reads
    rows = cursor.fetchall()
    for line in rows:
        total_count+=1
        # @rows contains only one value
        consensus_id = line[0]
# Prepare statement: get all the reads for this consensus_id
        sel_reads_sth = ('SELECT seq_id, orient, seq \
      FROM {}.reads \
      WHERE consensus_id = "{}"'.format(conf["database"],consensus_id))
        # only build consensus, if more then n reads
        cursor.execute(sel_reads_sth)
        if cursor.rowcount >= int(n_consensus):
            cons_count += 1
        # create output fasta
            filename="cons_{}_seqs.cfasta".format(consensus_id)
            fullpath = os.path.join(output,filename)
            outfile=open(fullpath, "w")
            reads_output = cursor.fetchall()
   # write selected reads to fasta file for muscle
            for lines in reads_output:
                seq_id = lines[0]
                orient =lines[1]
                seq= lines[2]
                if orient== "R":                # if necessary, reverse complement
                    seq=revcomp(seq)
                        # print each read to the fasta file
                test = ">{}\n{}\n".format(seq_id,seq)
                outfile.write(test)
            outfile.close()
### 3. Log information on consensi counts
    print("In total" , total_count , "of consensus_ids without sequence where found." , "\n",cons_count, "of them had more than", n_consensus, "reads and where thus considered for multiple sequence alignment.", "\n")
