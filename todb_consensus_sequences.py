#!/usr/bin/python3
"""
@Srilakshmy - s.harikrishnan@dkfz-heidelberg.de
@franasa - f.arcila@dkfz-heidelberg.de
    -adapted from Katharina Imkeller-

DESCRIPTION

Read the alignment output of MUSCLE, write consensus and identity percentage for each position to database.

1. Prepare database: insert consensus to sequences, update seq_id in consensus_stats
2. Get Consensus from alignment. Consensus file must have the name matching 'cons_[consensus_id]'.
3. Find out whether it is 1. or 2. consensus rank
4. Write to database.

LOGGED INFORMATION todo

- sequence ambiguities with keys and values (if not more then 50% of sequences have the one nucleotide)
- SQL statement: select all consensi refering to that event
- insert sequence statement
- SQL statement: update sequence_id in consensus_stats table

# Created: August 26th 2019

"""

import argparse
from argparse import RawTextHelpFormatter
from Bio import AlignIO
import os
import bcelldb_init as binit
import collections

##usage -> todb_consensus_sequences.pl -aln <musclealign> [-h]
##get configuartion for scireptor

config_scireptor=binit.get_config()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-aln', '--alignment', required=True, help="define input file")
    args = parser.parse_args()
    align = args.alignment
### 1. Prepare database for insertion of sequence
    conn = binit.connect()
##create a cursor object
    cursor = conn.cursor()
###prepare insert command


    # find out, whether it is 1. or 2. consensus
    # select the row, col, locus and experiment_id and check, whether there is another one for that

    sel_consensus_details = "SELECT row_tag, col_tag, experiment_id, locus \
      FROM {}.consensus_stats WHERE consensus_id = %s;".format(config_scireptor['database'])

    sel_all_consensi = 'SELECT consensus_id, n_seq \
  FROM {}.consensus_stats \
  WHERE row_tag = "%s" AND col_tag = "%s" AND experiment_id = "%s" AND locus = "%s" \
  ORDER BY n_seq DESC;'.format(config_scireptor['database'])

    # update consensus table, put seq_id
    update_consensus_st = "UPDATE {}.consensus_stats \
      SET sequences_seq_id='%s' WHERE consensus_id='%s'; ".format(config_scireptor['database'])

### 2. get consensus from alignment
    aln_file=os.path.basename(align).split(".")[0]
    position_count_hash = {}
    name_mapping = {}
    (cons, cons_id)= (aln_file.split("_")[0], aln_file.split("_")[1])

    aln_str = AlignIO.read(open(align), "fasta") # changed from perl str bc overwrites the method str

    cons_length = aln_str.get_alignment_length()


    for record in aln_str:

        sequence=record.seq

    for i in range(0,cons_length):

        nt_count = collections.Counter(aln_str[:, i])
        position_count_hash[i]= nt_count

    # find the most prominent letter
    # if it is not A,C,G or T, dont write it
    consensus= []
    consensus_quality = []

    for i in range(cons_length):

        max_value= (position_count_hash[i].most_common(1)[0][1])
        max_value_perc = (100 * max_value /len(aln_str))
        nt_ = (position_count_hash[i].most_common(1)[0][0])

        if nt_ in ["G","A","T","C"] and max_value >= 0.5*len(aln_str):

            consensus.append(nt_)
            consensus_quality.append(str(int(max_value_perc)))

consensus=("".join(consensus))
consensus_quality=(" ".join(consensus_quality))

### 3. Find out whether it is 1. or 2. consensus rank

consensus_rank = 1

cursor.execute(sel_consensus_details, [cons_id])

row, col, experiment_id, locus = cursor.fetchall()[0]


cursor.execute(str(sel_all_consensi % (row, col, experiment_id, locus)))

# log statement
# print ("\nSelect all consensi statement: %s\n With values (%s,%s,%s,%s).\n\n" % (sel_all_consensi,
#                                                                                 row, col, experiment_id, locus))
poss_cons_id, poss_nseq = cursor.fetchall()[0]

if not int(poss_cons_id) == int(cons_id):
    consensus_rank+=1 # todo works but compare to .pl expected behaviour


ins_seq_st = "INSERT INTO {}.sequences (seq, length, quality, name, consensus_rank) \
VALUES ('{}',{},'{}','{}',{}) ON DUPLICATE KEY \
UPDATE seq='{}' AND length={} AND quality='{}' AND name={} AND consensus_rank={} AND seq_id = LAST_INSERT_ID(seq_id);" \
    .format(config_scireptor['database'],str(consensus),len( consensus), consensus_quality, cons_id,
                consensus_rank, str(consensus),len( consensus), consensus_quality, cons_id, consensus_rank) #todo back to %s


try:
    cursor.execute(ins_seq_st)
except binit.MySQLdb._exceptions.OperationalError:
    print("{} was already inserted".format(ins_seq_st))


seq_id=cursor.lastrowid



cursor.execute(update_consensus_st %  (seq_id, cons_id))