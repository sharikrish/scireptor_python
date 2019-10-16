#!/usr/bin/python

####################################
# File name: todb_CDR_FWR.py#
# Author: Srilakshmy               #
# Created Date: August 29th 2019   #
####################################
import argparse
import bcelldb_init as binit
import re
from Bio import Seq
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
##get configuartion for scireptor

config_scireptor=binit.get_config()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-io', '--igblast_out', required=True, help="define input file")
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    igblastout = args.igblast_out

###Prepare database for insertion of sequence
### Logging and database init
    conn = binit.connect()
##create a cursor object
    cursor = conn.cursor()

# insert into warnings
    ins_warnings = "INSERT IGNORE \
  INTO {}.warnings (seq_id, FWR3_igblast_output, CDR3_start_C, CDR3_end, alt_CDR3_end, J_end) \
  VALUES (?,?,?,?,?,?);".format(config_scireptor['database'])

### define subroutines
    seq_id_sequence = {}
    seq_id_locus={}

    def revcomp(dna):
        seq_dna=Seq(dna)
        return seq_dna.reverse_complement()

    def get_regions(type, start,end, seq_id):
        sequence = seq_id_sequence[seq_id]
        if(end-start+1 > 0):
            dna_seq= sequence[start-1:(end-start+2)+start]
            nucobj = Seq(dna_seq,generic_dna)
            protseq= nucobj.translate()
            prote_seq = str(protseq)
            protlen = len(protseq)
            if prote_seq.find("*") != -1:
                stop_codon=1
            else:
                stop_codon=0

            # insert into CDR_FWR table
            ins_regions = "INSERT IGNORE \
          INTO {}.CDR_FWR (seq_id, region, start, end, dna_seq, prot_seq, prot_length, stop_codon) \
          VALUES (%s,%s,%s,%s,%s,%s,%s,%s);".format(config_scireptor['database'])
            val= (seq_id, type, start, end, dna_seq, prote_seq, protlen, stop_codon)
            cursor.execute(ins_regions, val)

        # second, assign the RegExp to the respective hash, using the locus as key. ATTENTION: key is upper case
    CDR3_end_motif ={}
    alt_CDR3_end_motif={}
    J_end_motif={}
    locus = ('a','b','h','k','l')
    for key_conf_locus in locus:
        # first, remove the quotes flanking the RegExp in the config file
        for key_conf_suffix in ('CDR3_e', 'altCDR3_e1', 'altCDR3_e2', 'altCDR3_e3', 'Jend' ):
            temp_key_conf=key_conf_locus + "_" + key_conf_suffix
            try:
                if config_scireptor[temp_key_conf] and config_scireptor[temp_key_conf] != '^""$':
                    config_scireptor[temp_key_conf]= re.search(r'(^.*)$',config_scireptor[temp_key_conf]).group(0)
                else:
                    config_scireptor[temp_key_conf] = "ZZZZ"
            except NameError:
                print("No")
        # second, assign the RegExp to the respective hash, using the locus as key. ATTENTION: key is upper case
        uc = key_conf_locus.upper()
        CDR3=key_conf_locus + "_CDR3_e"
        CDR3_end_motif[uc] = config_scireptor[key_conf_locus + "_CDR3_e"]
        alt_CDR3_end_motif[uc]=(config_scireptor[key_conf_locus + "_altCDR3_e1"],
                                config_scireptor[key_conf_locus + "_altCDR3_e2"],
                                config_scireptor[key_conf_locus + "_altCDR3_e3"])
        J_end_motif[uc]=config_scireptor[key_conf_locus + "_Jend"]
#
#     # Generate hashes containing the CDR3 and J motifs. Each locus has its own set of motifs and up to three alternative
#     # CDR3 motifs are supported. As noted in the config file, the motif must exist and must not be empty, otherwise the
#     # resulting RegExp ("//") would simply match the first character. This is assumed to be an unindented behavior, thus
#     # in these cases the key will be assigned the irrelevancompare(a1,a2,allowAll=TRUE)t pattern "ZZZZ", which does not exist in amino acid strings.
#     ###
#     ### 2. get query_ids and CDR/FWR positions from IgBLAST output
    count_line=0
    count_query=0
    summary_mark=0
    with open(igblastout) as ig:
        lines=ig.readlines()
        for i, line in enumerate(lines):
            count_line+=1
###############################################################################
            if "Query:" in line:
                count_query+=1
                # split the query line to extract the id
                (comment, title, query_id) = line.split(" ")
                query_id=query_id.rstrip("\n")
                # reset the variables
                # needed, in case IgBLAST output not complete for certain query
                (type, start, end)= ["N/A"] * 3
                (FWR3_start, FWR3_end, CDR3_start, CDR3_end)= ["N/A"] * 4
                # get sequence from database
                sel_sequence = 'SELECT locus, orient, seq \
                  FROM {}.sequences WHERE seq_id = "{}";'.format(config_scireptor['database'], query_id)
                cursor.execute(sel_sequence)
                seq_array=cursor.fetchall()
                for row in seq_array:
                    seq_id_locus[query_id]= row[0]
                    orient = row[1]
                    if(orient == "R"):
                        seq_id_sequence[query_id] = revcomp(row[2])
                    elif(orient== "F"):
                        seq_id_sequence[query_id] = row[2]
                    query_length = len(seq_id_sequence[query_id])
###############################################################################
            elif "Alignment summary" in line:
                summary_mark=count_line
##################################################################################
            ### 3. extract the positions from the alignment summary
            elif count_line >= summary_mark+1 and count_line <= summary_mark+6 and count_line >10:
                all_list = line.split("\t")
                if len(all_list) >1:
                    (type, start, end) = (all_list[0], all_list[1], all_list[2])
                    type= re.search(r'^(CD|F)R[1-3]', type)
                    if type:
                        type= type.group(0)
            # write the positions and sequences for regions until CDR2 into output
                simple_regions = ("FR1", "CDR1", "FR2", "CDR2")
                if str(type) in line and str(type) in simple_regions:
                        if not start == "N/A":
                           get_regions(type,int(start),int(end),query_id)
                if str(type) =="FR3":
            # store the positions of FWR3
            # since looking for CDR3 is more difficult, it will be done separately
            # and using the information about FWR3
                    FWR3_start=start
                    FWR3_end = end
                ### 4. Split sequence in regions and write to database
                # At the begin of a new query block process the remainder of the previous query sequence to extract positions
                # of FWR3, CDR3, J and constant. (count_line >= 3) is required to ignore the first lines of the IgBLAST output,
                # "# BLAST processed ..." to process the last query block (since "# IGBLASTN ..." starts a new block and is
                # therefore absent for the last one.
                # CDR3 positions have not been initialized yet, since they do not appear in the IgBLAST output
                # or we like to define them differently
                # initialize variables for warnings table
                elif r"# IGBLASTN [0-9].+$" and count_line >3 or r"# BLAST processed [0-9]+ quer(?:y|ies)":
                    FWR3_igblast_output = 0
                    CDR3_start_C = 0
                    CDR3_end = 0
                    alt_CDR3_end= 0
                    rest_J_protstart = 0
                    # initialize variables for warnings table
                    if not (FWR3_start == "N/A"):
                        FWR3_igblast_output=1
                        # find region that should contain the CDR3
                        FWR3dna_length=int(FWR3_end)-int(FWR3_start) + 1
                        frame_overhang = FWR3dna_length % 3
                        critregion_start = int(FWR3_start) + (FWR3dna_length - frame_overhang - 15)
                        sequence= seq_id_sequence[query_id]
                        critregion=sequence[critregion_start-1:critregion_start+155]
                        nucleo = Seq(critregion, generic_dna)
                        protseq = nucleo.translate()
                        critprot = str(protseq)

                # find actual CDR3 in terms of aa
                        (CDR3protend,CDR3protstart) = ("N/A", "N/A")
                        if "C" in critprot:
                            CDR3protstart= re.search("C", critprot).end()
                            CDR3_start_C=1
                        # go on with looking for CDR3 end motif
                        (rest_J_protstart, rest_J_protend, rest_const_protstart, rest_const_protend) = ("N/A", "N/A", "N/A", "N/A")
                        regex=re.compile(CDR3_end_motif[seq_id_locus[query_id]])
                        if regex.search(critprot):
                            CDR3protend = regex.search(critprot).end()
                            CDR3protend=CDR3protend-5
                            rest_J_protstart = CDR3protend + 1
                            CDR3_end=1
                        for m in alt_CDR3_end_motif[seq_id_locus[query_id]]:
                            regex=re.compile(m)
                            for match in regex.finditer(critprot):
                                CDR3protend= match.end()
                                CDR3protend = CDR3protend - 5
                                rest_J_protstart = CDR3protend + 1
                                alt_CDR3_end = 1
                # look for necessary AA motif and adjust rest protein positions
                        regex = re.compile(J_end_motif[seq_id_locus[query_id]])
                        if regex.search(critprot):
                            rest_const_protstart = regex.search(critprot).end()
                            rest_J_protend = rest_const_protstart - 1
                            J_end = 1
                        if CDR3_end_motif == 0 and (alt_CDR3_end == 0 or J_end == 0):
                            CDR3protend = "N/A"
                            rest_J_protstart = "N/A"
                        if(CDR3protstart != "N/A" and CDR3protend!="N/A"):
                            CDR3_end_pos = CDR3protend-CDR3protstart+1
                            CDR3_end_pos=CDR3protstart+CDR3_end_pos
                            CDR3=critprot[CDR3protstart:CDR3_end_pos]
                            CDR3_length= len(CDR3)
                            CDR3dna_length=CDR3_length*3
                            # reset FWR3 end if necessary
                            FWR3_end = critregion_start + 3 * CDR3protstart - 1
                            CDR3_start = FWR3_end + 1
                            CDR3_end = CDR3_start + CDR3dna_length - 1
                            get_regions("CDR3", int(CDR3_start), int(CDR3_end), query_id)

                        # write regions for FWR3 into output
                        # this is done here, due to optional reseting, after CDR3 was localized

                        get_regions("FR3", int(FWR3_start), int(FWR3_end), query_id)

                        if(rest_J_protstart!="N/A" and rest_J_protend!= "N/A"):
                            rest_J_start = critregion_start + 3 *rest_J_protstart
                            rest_J_end = critregion_start + 3 * (rest_J_protend+1) - 1
                            get_regions("rest_J", int(rest_J_start), int(rest_J_end), query_id)

                        if (rest_const_protstart!= "N/A" and rest_J_protend != "N/A"):
                            rest_const_start = critregion_start + 3 *rest_const_protstart
                            rest_const_end = query_length
                            get_regions("rest_const", int(rest_const_start), int(rest_const_end), query_id)




















