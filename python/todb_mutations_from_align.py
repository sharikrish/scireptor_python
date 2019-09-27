#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de  -adapted from Katharina Imkeller-

DESCRIPTION

This function calculates the mutations from the query-germline alignment and writes them to the database.
Take all alignments in input_directory with information on start of query and germline in reference to origin (frame 0). Take the mutation matrix and check for every codon (in frame!), where replacements, silent mutations, ins and dels occur. The variables in_status and del_status tracks the relative frame shift in reference to origin. BE CAREFULL: If the frame is shifted due to insertion/deletion, the translation of single codons is not meaningfull anymore and you have to check, whether silent/replacement assignment really correspond to what you want to see!
In the alignment file it allways needs to be like that:

><queryid>_<querystartposition>_query
<query sequence>
><queryid>_<germlinestartposition>_germline
<germline sequence>

1. Get all files from input directory
2. Build up hash from mutation_matrix
3. Go through all files, split alignment into reading frame codons
4. For every codon look for mutations and stopcodons

TODO LOGGED INFORMATION

- sequences where there is a stop codon either in the query or in the germline

TODO

deal with letters different than ACGT

optionaly write the output into a file (if just needed for analysis and do not want to write to DB)

"""


import argparse
from argparse import RawTextHelpFormatter
import os
import bcelldb_init as binit
from Bio import AlignIO, SeqIO

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-dir', '--directory', required=True, help="define output file")
args = parser.parse_args()

inputdir = args.directory


conf = binit.get_config()

db = binit.connect()
cursor = db.cursor()


## logging to database currently turned off due to database time-outs

###
### 1. Get all files from input directory
###
# use all input files ending on .igblast.aln in the directory

aln_files = [os.path.join(inputdir,aln) for aln in os.listdir(inputdir) if aln.endswith(".igblast.aln")]

mutation_table_hash = open("mutation_matrix.txt","r")


mutation_dict={}

for line in mutation_table_hash.readlines():
    obs_codon,germ_codon,mutation,n_repl,n_silent,insertion,deletion = (line.strip().split("\t"))
    # mutations without known effect (that additionally occur in a codon with in del)
    mutation_dict[(obs_codon, germ_codon, "mutation")] = mutation
    # replacement mutations
    mutation_dict[(obs_codon, germ_codon, "n_repl")] = n_repl
    # silent mutations
    mutation_dict[(obs_codon, germ_codon, "n_silent")] = n_silent
    # insertions
    mutation_dict[(obs_codon, germ_codon, "insertion")] = insertion
    # deletions
    mutation_dict[(obs_codon, germ_codon, "deletion")] = deletion



# print (aln_files)




###
### 3. Go through all files, split alignment into reading frame codons
###


for file in aln_files:

    records = list(SeqIO.parse(file, "fasta"))
    # query: nr. 1, subj: nr. 2
    query_seq = records[0].seq

    # determine, where on the alignment to start with the reading frame
    # split headers of query and germline
    # extract starting positions
    subject_seq = records[1].seq
    query_start = int((records[0].id).split("_")[1])
    germline_start = int((records[1].id).split("_")[1])

    # determine frame and find out, where first codon starts
    frame = germline_start % 3
    # the "frame" does not correspond to the position, where the first codon starts
    # need the variable $frame_start to
    if frame == 2:
        frame_start = 2
    elif frame == 1:
        frame_start = 0
    elif frame == 0:
        frame_start = 1

    # consecutive_in_del_status
    # +1 for insertions
    # -1 for deletions
    in_status = 0
    del_status = 0

    statement = 'INSERT IGNORE INTO {}.mutations ' \
                '(seq_id, position_codonstart_on_seq, replacement, silent,' \
                'insertion, deletion, undef_add_mutation,' \
                'stop_codon_germline, stop_codon_sequence,' \
                'in_status, del_status) ' \
                'VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'.format(conf["database"])

    seq_id = int(file.split("/")[-1].split('.igblast.aln')[0])



    # print seq_id into log table
    if int(conf["log_level"]) >= 4:
        print(f"[todb_mutations_from_align.pl][DEBUG] Processing sequence {seq_id}\n")

    # print(frame, seq_id)

    ###
    ### 4. For every codon look for mutations and stopcodons
    ###
    # split the alignments into codons and check for mutations
    # start at the position previously determined

    for i in range(frame_start, len(str(query_seq)), 3): # todo change to alligment lenght

        query_codon = str(query_seq)[i:i + 3]
        subject_codon = str(subject_seq)[i:i + 3]

        # print(subject_codon)

        # reset boolean variables to check for stop codons
        stop_codon_germline = 0
        stop_codon_sequence = 0

        # reset mutation variables
        replacements = 0
        silents = 0
        add_mutations = 0
        insertions = 0
        deletions = 0

        # is there a stop codon in the germline? --> then there likely is a problem with the frame...
        if subject_codon == 'TAA' or subject_codon == 'TGA' or subject_codon == 'TAG':
            if int(conf["log_level"]) >= 4:
                print (f'[todb_mutations_from_align.pl][INFO] sequence $seq_id position {i} query \ '
                       f'{query_codon} \ germline \ {subject_codon} \  : stop codon in germline sequence.\n')
            stop_codon_germline = 1

        if query_codon == 'TAA' or query_codon == 'TGA' or query_codon == 'TAG':
            if int(conf["log_level"]) >= 4:
                print(f'[todb_mutations_from_align.pl][INFO] sequence $seq_id position {i} query \ '
                      f'{query_codon} \ germline \ {subject_codon} \  : stop codon in query sequence.\n')
            stop_codon_sequence = 1

        # try to find mutations only, if the two codons are different, # make sure the codon is 3 nucl long
        # also insert when there is any stop codon
        if subject_codon != query_codon and len(query_codon) == 3 or stop_codon_sequence == 1 or stop_codon_germline == 1:
            # print (query_codon, subject_codon)
            if mutation_dict[(query_codon, subject_codon, "mutation")]:
                print(mutation_dict[(query_codon, subject_codon, "mutation")], subject_codon, query_codon)
                replacements = mutation_dict[(query_codon, subject_codon, "n_repl")]
                silents = mutation_dict[(query_codon, subject_codon, "n_silent")]
                add_mutations = mutation_dict[(query_codon, subject_codon, "mutation")]
                insertions = mutation_dict[(query_codon, subject_codon, "insertion")]
                deletions = mutation_dict[(query_codon, subject_codon, "deletion")]

                in_status += int(insertions)
                del_status += int(deletions)

            # determine the position of the mutation in terms of absolute measure on the query sequence

            position_codon_on_seq = i + query_start

            # insert into mutations table of database
            val = (seq_id, position_codon_on_seq, replacements, silents, insertions, deletions,
                   add_mutations, stop_codon_germline, stop_codon_sequence, in_status, del_status)

            temp_insert_return = cursor.execute(statement % val)

            if not temp_insert_return:
                if int(conf["log_level"]) >= 1:
                    print(f"[todb_mutations_from_align.pl][ERROR] Could not insert into table. {statement % val}->errstr \n")
            else:
                if int(conf["log_level"]) >= 4:
                    print(f"[todb_mutations_from_align.pl][DEBUG] Inserted {temp_insert_return} lines into table.\n")




            # print (position_codon_on_seq)





mutation_table_hash.close()