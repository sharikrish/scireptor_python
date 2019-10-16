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
from math import log

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
    #hash_ascii_phred = {"!": "0", '"' : "1", "#": "2", "$": "3", "%": "4", "&": "5","'": "6","(": "7",")":"8", "*": "9","+": "10",",":"11","-": "12", "."	: "13","/":"14", "0":"15","1": "16","2": "17","3":  "18", "4": "19", "5": "20", "6": "21","7":"22", "8": "23", "9":"24",":":"25",";": "26","<": "27","=": "28", ">"	: "29", "?" : "30", "@"	: "31", "A"	: "32", "B"	: "33", "C"	:"34", "D"	: "35", "E"	: "36", "F"	: "37", "G"	: "38", "H"	: "39", "I": "40"}
    hash_ascii_phred = {'0': '!', '1': '"', '2': '#', '3': '$', '4': '%', '5': '&', '6': "'", '7': '(', '8': ')', '9': '*', '10': '+', '11': ',', '12': '-', '13': '.', '14': '/', '15': '0', '16': '1', '17': '2', '18': '3', '19': '4', '20': '5', '21': '6', '22': '7', '23': '8', '24': '9', '25': ':', '26': ';', '27': '<', '28': '=', '29': '>', '30': '?', '31': '@', '32': 'A', '33': 'B', '34': 'C', '35': 'D', '36': 'E', '37': 'F', '38': 'G', '39': 'H', '40': 'I'}

    sel_consensus_details = 'SELECT row_tag, col_tag, experiment_id, locus \
      FROM {}.consensus_stats WHERE consensus_id = %s;'.format(config_scireptor['database'])

    sel_all_consensi = 'SELECT consensus_id, n_seq \
      FROM {}.consensus_stats \
      WHERE row_tag = "%s" AND col_tag = "%s" AND experiment_id = "%s" AND locus = "%s" \
      ORDER BY n_seq DESC;'.format(config_scireptor['database'])

    # update consensus table, put seq_id
    update_consensus_st = "UPDATE {}.consensus_stats \
        SET sequences_seq_id='%s' WHERE consensus_id='%s'; ".format(config_scireptor['database'])

### 2. get consensus from alignment

    aln_file = os.path.basename(align).split(".")[0]
    position_count_hash = {}
    name_mapping = {}
    cons_id = aln_file.split("_")[1]
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
    (old_min, old_max,new_min,new_max)= (50,100,0,40)
    for i in range(cons_length):
        max_value = (position_count_hash[i].most_common(1)[0][1])
        nt_ = (position_count_hash[i].most_common(1)[0][0])
        if nt_ in ["G","A","T","C"] and max_value >= 0.5*len(aln_str):
            max_value_perc = (100 * max_value / len(aln_str))
            Qscore = int((-10 * log(10 ** (-max_value_perc / 10), 10)))
            new_value = (int((((Qscore - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min))
            if new_value < 0:
                print(cons_id, "\n", max_value, new_value, position_count_hash[i].most_common(1)[0][1])
            ascii_value = hash_ascii_phred[str(new_value)]
            consensus.append(nt_)
            consensus_quality.append(str(ascii_value))
consensus=("".join(consensus))
consensus_quality=(("".join(consensus_quality)))
### 3. Find out whether it is 1. or 2. consensus rank
consensus_rank = 1
cursor.execute(sel_consensus_details, [cons_id])
row, col, experiment_id, locus = cursor.fetchall()[0]
cursor.execute(str(sel_all_consensi % (row, col, experiment_id, locus)))
# log statement

poss_cons_id, poss_nseq = cursor.fetchall()[0]
if not int(poss_cons_id) == int(cons_id):
    consensus_rank+=1 # todo works but compare to .pl expected behaviour

# a combination of placeholder %s and formatting {} is used to scape quotes and single quotes from quality strings

ins_seq_st = "INSERT INTO {}.sequences (seq, length, quality, name, consensus_rank) \
VALUES ('{}',{},%s,'{}',{}) ON DUPLICATE KEY \
UPDATE seq='{}' AND length={} AND quality=%s AND name={} AND consensus_rank={} AND seq_id = LAST_INSERT_ID(seq_id);" \
    .format(config_scireptor['database'],str(consensus),len( consensus), cons_id,
                consensus_rank, str(consensus),len( consensus), cons_id, consensus_rank)

try:
    cursor.execute(ins_seq_st,(consensus_quality, consensus_quality))
except binit.MySQLdb._exceptions.OperationalError:
    print("{} was already inserted".format(ins_seq_st))
# except binit.MySQLdb._exceptions.ProgrammingError:
#     print (cons_id, ins_seq_st, binit.MySQLdb._exceptions.ProgrammingError)



seq_id=cursor.lastrowid

cursor.execute(update_consensus_st % (seq_id, cons_id))