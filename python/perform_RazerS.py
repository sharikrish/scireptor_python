#!/usr/bin/python

####################################
# File name: perform_RazerS.py     #
# Author: Srilakshmy               #
# Created Date: August 15th 2019     #
####################################

import argparse, subprocess
import bcelldb_init as binit
import os

##get configuartion for scireptor

config_scireptor=binit.get_config()
conn=binit.connect()

##create a cursor object
cursor = conn.cursor()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="define input file")
    parser.add_argument('-ro', '--output', required=True, help="define output file")
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

###define input and output file	###

    filename = args.fasta
    output = args.output
    dir, extension = os.path.split(output)
    tag_name = extension.split('.')[0]

    ###### 1. Get tags from library and write to fasta file
    # get tags from the database

    get_tag_statement = 'SELECT tag_id, sequence \
      FROM {}.tags_library \
      WHERE matrix=\"{}\" AND batch=\"{}\";' \
        .format(config_scireptor['library'], config_scireptor['matrix'], config_scireptor['tag_batch'])

    # log select
    print('\n Tag select statement: ' + get_tag_statement + '\n\n')

    number_of_rows = cursor.execute(get_tag_statement)
    result = cursor.fetchall()
    input_tag =os.path.join(dir, str(tag_name) + '.seqfasta.tags.fa')
    tag_fasta_outfile = open(input_tag, 'w')
    for row in result:
        tag_fasta_outfile.write(">" + str(row[0]) + "\n" + str(row[1]) + "\n")
    tag_fasta_outfile.close()
    ### 2. Perform razers

# razers parameters
    razers_percid="90"    # minimal required percid between query and template
    m="100000000"        # maximal number of displayed results
    s="1111111"          # shape factor, seven ones mininal, means that there need to be seven in a row
    pf="1"          # way to present positions in string, 1 => first position=1 not 0
###call Razer###
    try:
        string_process= subprocess.call(["/usr/bin/razers3", "-a", "-m", m, "-s", s, "-i", razers_percid, "-o", output, "-pf", pf, filename, input_tag])
    except OSError:
        print('\nCould not find Razers3, Please reinstall. \n')


