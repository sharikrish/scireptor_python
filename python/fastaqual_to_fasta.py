#!/usr/bin/python3

"""
Convert FASTA + QUAL file pairs to a single FASTQ file
http://seqanswers.com/forums/showthread.php?t=16925

You can use this script from the shell like this::
$ ./fasta_to_fastaq reads.fna reads.qual reads.fastq
"""

# The libraries we need #
import sys, os,argparse
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
# Get the shell arguments #
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True, help="define input file")
    parser.add_argument('-q', '--qual', required=True, help="define quality file")
    parser.add_argument('-p', '--output_path', required=True, help="define output file")
    args = parser.parse_args()
    fa_path=args.fasta
    qa_path=args.qual
    fq_path = args.output_path

    # create output fasta
    base=os.path.basename(fa_path)
    print(base)
    filename = base + ".fastq"
    fullpath = os.path.join(fq_path, filename)
    outfile = open(fullpath, "w")
# # Do it #
    records = PairedFastaQualIterator(open(fa_path), open(qa_path))
    count = SeqIO.write(records, outfile, "fastq")
    outfile.close()
# # Report success #
print("Converted %i records" % count)
