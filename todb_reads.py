#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de  -adapted from Katharina Imkeller-

DESCRIPTION 

Load the raw high-throughput sequences into the database, including quality 
information and information about sequencing run.

1.	Get sequences and qualities from input files and store them into hashes with identifier as keys.
    Create a list of all identifiers.
2.	Prepare statements to insert into reads and sequencing_run tables.
3.	Write sequencing run info into sequencing run table and get back the identifier.
4.	For each identifier write sequence, length, quality and sequencing_run identifier into reads table.

LOGGED INFORMATION

- information on insertion of sequencing_run
- number of reads that were inserted
"""


import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import bcelldb_init as binit
import re
from math import log

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument("-fq", "--fastq", required=True,  help="Fasta file")
#parser.add_argument("-q","--qualfile", required=True,  help="Quality file")
parser.add_argument("-ri","--seqrun_info", required=True,  help="Text file with information concerning the sequencing run")
parser.add_argument("-m","--matrix", required=True,  help="Experiment ID (formerly called matrix)")
args = parser.parse_args()

config_seq_length_max = 1000
# output = "reads"  # not used variable
fastafile = args.fastq
#qualfile = args.qualfile
seqrun_info = args.seqrun_info
matrix = args.matrix

### 0. Logging information

db = binit.connect()
cursor = db.cursor()
print(fastafile)


### 1. Get sequences and qualities from input files
### and store them into hashes with identifier as keys.
### Create a list of all identifiers.
### Does not make sure that for every identifier there exist sequence and quality!


# open file with sequence information

fasta_in =SeqIO.index(fastafile, "fastq")
#qual_in = SeqIO.index(qualfile, "qual")

conf = binit.get_config()

# identifier list
identifiers=list(fasta_in.keys())

### 2. Get sequencing run from infile and put it into a dictionary


seqrun_infofile = open(seqrun_info, "r")

infile_dic={}
for line in seqrun_infofile:
    if not line.startswith("#") and not re.match(r'^\s*$', line):
        key= line.strip().split("=")[0]
        val = line.strip().split("=")[1]
        infile_dic[key] = val


### 3. Prepare statements to insert into reads and sequencing_run tables
### and write into sequencing run table
### Get back the identifier.




statement2 = 'INSERT IGNORE INTO {}.sequencing_run(' \
             'date, name, processed_by, plate_layout_id,' \
             ' add_sequencing_info, experiment_id)VALUES (' \
             '"{}","{}","{}","{}","{}","{}")'.format(conf['database'],infile_dic['rundate'],  #
                                         infile_dic['runname'],
            infile_dic['processed_by'], infile_dic['platelayout_id'],
                                         infile_dic['optional'], matrix)
# execute query


print (infile_dic['rundate'])

seq_run_bool = cursor.execute(statement2)
sequencing_run_id = db.insert_id()


# print(seq_run_bool)

# db.insert_id()

# get the last inserted id, which corresponds to sequencing_run_id


if seq_run_bool == True:
    sequencing_run_id = db.insert_id()
    sequencing_run_id = int(cursor.lastrowid) # check this behaviour later
    print ("[todb_reads.py][INFO] Inserted run {} {} from {} as id {} into table \"sequencing_run\".\n".format(
        infile_dic['runname'], infile_dic['optional'], infile_dic['rundate'], "sequencing_run_id"
    ))
    pass
else:
    print ("[todb_reads.py][WARNING] Run {} from {} was already present in table {} and not inserted.\n".format(
        infile_dic['runname'],infile_dic['optional'], infile_dic['rundate']
    ))



### 4. For each identifier write sequence, length,
### quality and sequencing_run_id into reads table.

# reset counter for the number of total and inserted reads
n_reads = 0
n_inserted = 0
def ave_qual(quals):
    """Calculate average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        return -10 * log(sum([10**(q / -10) for q in quals]) / len(quals), 10)
    else:
        return None

for id in identifiers:
    n_reads+=1
    id_seq = fasta_in[id].id[0:config_seq_length_max]
    seq = fasta_in[id].seq  # gets the raw sequence
    qual_id= fasta_in[id].format("fastq").split('+')[1].strip('\n')

    #list_qual_id = (qual_in[id].letter_annotations["phred_quality"])[0:config_seq_length_max]
    #qual_id = (" ".join(str(e) for e in list_qual_id))

    statement1 = 'INSERT IGNORE INTO {}.reads (' \
                 'name, length, seq, quality, ' \
                 'sequencing_run_id) VALUES (' \
                 '"{}","{}","{}","{}","{}") '.format(conf["database"], id, len(seq),
                                                     seq, qual_id, sequencing_run_id)

    cursor.execute(statement1)

print ("total number of reads = {}".format(str(n_reads)))


# TBD: LOG check for raise errors, check what happens when running a second time
