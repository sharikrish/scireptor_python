#!/usr/bin/python

####################################
# File name: perform_muscle.py     #
# Author: Srilakshmy               #                 
#Created Date: August 8th 2019     #
####################################

import sys,argparse, subprocess

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta', required=True, help="define input file")
	parser.add_argument('-aln', '--output', required=True, help="define output file")
	parser.add_argument('-v', dest='verbose', action='store_true')
	args = parser.parse_args()
        filename = args.fasta
        output = args.output
        
###call muscle###
        try:
            subprocess.call(["/usr/bin/muscle", "-in", filename, "-maxiters", "1" , "-quiet", "-out", output])
        except OSError:
            print('\nCould not find muscle, Please reinstall. \n')
