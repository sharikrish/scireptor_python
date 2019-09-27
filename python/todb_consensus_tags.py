#!/usr/bin/python

####################################
# File name: todb_consensus_tags.py#
# Author: Srilakshmy               #
# Created Date: August 16th 2019   #
####################################
import dbconnect
import sys, argparse, subprocess
import bcelldb_init as binit
from configobj import ConfigObj
import re
seqid_orient_hash= {}
seqid_locus_hash= {}
wellid_hash={}
seq_list=[]
hash_locus_num = {"H" : "1",  "K" : "2",  "L"  : "3",  "B"  :"4",  "A" : "5"}
hash_num_locus = {'1' : "H", '2' :"K", '3' : "L", '4' : "B", '5' : "A"}

##get configuartion for scireptor

config_scireptor=binit.get_config()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--experiment_id', required=True, help="define input file")
    parser.add_argument('-l', '--locus', required=True, help="define output file")
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    experiment_id = args.experiment_id
    fixed_locus = args.locus

###Prepare database for insertion of sequence
### Logging and database init
    conn = dbconnect.connect()
##create a cursor object
    cursor = conn.cursor()

    ### 1. Get a list of all seq_id that will be processed,
    ### i.e. all the ones that do not yet have a well id assigned.
    # the ones that do not yet have a consensus but at least two tags

    sel_reads1_sth= "SELECT reads.seq_id, orient, locus FROM %s.reads\
	JOIN \
	(SELECT reads_tags.seq_id, COUNT(reads_tags.seq_id) AS cnt FROM %s.reads_tags \
	GROUP BY reads_tags.seq_id) AS count_tags\
	ON count_tags.seq_id = reads.seq_id\
	WHERE reads.well_id IS NULL\
	AND locus = '%s'\
	AND orient IS NOT NULL\
	AND count_tags.cnt >=2;" \
        %(config_scireptor['database'],config_scireptor['database'],fixed_locus)
    print(sel_reads1_sth)
    cursor.execute(sel_reads1_sth)
    result = cursor.fetchall()
    count_seqs=0
    for seq in result:
        seqid_orient_hash[seq[0]] = seq[1]
        seqid_locus_hash[seq[0]] = seq[2]
        seq_list.append(seq[0])
        count_seqs+=1
        # log how many sequences were found
    print "Selected", count_seqs, "reads for processing from locus", fixed_locus,"\n"
# ### 2. Prepare database query: for a certain seq_id get the right identifying tags from the database
# # select all the tags from the database
# # forward tag
# # looking for max(percid) and min(start)
    hash_wellid_seqid= {}
    for seqid in seq_list:
            sel_Ftags_sth = "SELECT found.tag_id, lib.name\
                FROM %s.reads_tags as found\
                JOIN %s.tags_library as lib\
                ON found.tag_id = lib.tag_id\
                WHERE seq_id ='%s' \
                AND direction = 'F'\
                ORDER BY percid DESC, start ASC\
                LIMIT 1; " % (config_scireptor['database'], config_scireptor['library'], seqid)

            sel_Rtags_sth = "SELECT found.tag_id, lib.name\
                FROM %s.reads_tags as found\
                JOIN %s.tags_library as lib\
                ON found.tag_id = lib.tag_id\
                WHERE seq_id ='%s' \
                AND direction = 'R'\
                AND start > '%s' \
                ORDER BY percid DESC, start DESC\
                LIMIT 1; " % (
            config_scireptor['database'], config_scireptor['library'], seqid, config_scireptor['tag_landing_zone'])
            cursor.execute(sel_Ftags_sth)
            ftag = cursor.fetchall()
            cursor.execute(sel_Rtags_sth)
            rtag = cursor.fetchall()
            if not (not rtag) or (not ftag):
            #     pass
            # else:
                for i in ftag:
                    (Fid, Ftagname) = (i[0], i[1])
                for j in rtag:
                    (Rid, Rtagname) = (j[0], j[1])
                if ((re.match("(^R)",Ftagname)) and (re.match("(^C)", Rtagname))):
                    # if forward is R and reverse is C, read orientation should be reverse
                    if (seqid_orient_hash[seqid] == "F"):
                        # log: tags and orientation do not match
                        print("Read", seqid,  "has a Ftag", Ftagname, "and a Rtag", Rtagname, "but an orientation ", seqid_orient_hash[seqid], "\n")
                    elif seqid_orient_hash[seqid] == "R":
                        # tags and orientation fit
                        # take only the number
                        (rowtag, coltag) = (Ftagname[1::], Rtagname[1::])
                    else:
                        print ("\nread",  seqid,  "has no correct orientation\n")	# log: no assigned orientation
                elif((re.match("(^C)",Ftagname)) and (re.match("(^R)", Rtagname))):
                # if forw is C and rev is R, read orientation should be forward
                    if (seqid_orient_hash[seqid] == "R"):
                        print("\nread", seq_d, "has a Ftag", Ftagname, "and a Rtag", Rtagname, ",but an orientation",  seqid_orient_hash[seq_id],  "\n")
                    elif seqid_orient_hash[seqid] == "F":
                        (rowtag, coltag) = (Rtagname[1::],Ftagname[1::])
                    else:
                        print("\nread 'seq_id' has no correct orientation\n")
                if (rowtag!= '' and coltag!= '') :
                    #print(rowtag)
                    wellid_hash = [coltag,rowtag, hash_locus_num[seqid_locus_hash[seqid]]]
                    hash_wellid_seqid[seqid]= wellid_hash
                    #print "read", seqid, "was matched to well_id", wellid_hash,  "\n";
    cnt_wellids_assigned = 0
    cnt_reads_updated = 0
    separtor =""
    for wellid_curr in hash_wellid_seqid.keys():
        seqid_curr = hash_wellid_seqid[wellid_curr]
        cnt_wellids_assigned+=1
        #cnt_reads_updated+= seqid_curr
        well_id_each=separtor.join(seqid_curr)
        statement_update_reads_wellid = "UPDATE %s.reads SET well_id = '%s' WHERE seq_id = '%s';" %(config_scireptor['database'], well_id_each, wellid_curr)
        cursor.execute(statement_update_reads_wellid)
    # ### 5. Now, select all occuring well_ids without consensus_id
    sel_wellid = "SELECT cast(well_id as char) FROM %s.reads\
                    WHERE consensus_id IS NULL\
                     AND well_id IS NOT NULL\
                     AND locus = '%s'\
                     GROUP BY well_id;" % (config_scireptor['database'],fixed_locus)
    cursor.execute(sel_wellid)
    ##adding counters
    count_wellids = 0
    count_seq_ids = 0
    wellid_without_con_id = cursor.fetchall()
    for k in wellid_without_con_id:
        count_wellids += 1
        well_id=k[0]
    # # convert well_id back to locus, row and col
        locus = hash_num_locus[well_id[-1]]
        rowtag = "R"+well_id[0:3]
        coltag = "C"+well_id[3:6]
        # 6. Prepare database: for a certain well_id, select the most common and 2nd most common V-J combination
        sel_VJs = "SELECT Vseg.VDJ_id, Jseg.VDJ_id, COUNT(*) as cnt \
                  FROM %s.reads \
                     JOIN %s.reads_VDJ_segments as Vseg ON reads.seq_id = Vseg.seq_id \
                     JOIN %s.reads_VDJ_segments as Jseg ON reads.seq_id = Jseg.seq_id \
                     JOIN %s.sequencing_run ON reads.sequencing_run_id = sequencing_run.sequencing_run_id \
                     WHERE Vseg.type='V' AND Jseg.type='J' \
                     AND Vseg.igblast_rank=1 AND Jseg.igblast_rank=1 \
                     AND reads.well_id = '%s' \
                     AND sequencing_run.experiment_id = '%s' \
                     GROUP BY Vseg.VDJ_id, Jseg.VDJ_id \
                     ORDER BY cnt DESC \
                     LIMIT 2;" %(config_scireptor['database'],config_scireptor['database'],config_scireptor['database'],config_scireptor['database'],well_id,experiment_id)
        # get most occuring V_J combinations
        cursor.execute(sel_VJs)
        best_VJ=cursor.fetchall()
    # ##7. update consensus_id of sequence with right well_id, V and J segment
        for vj in best_VJ:
            vseg=vj[0]
            jseg=vj[1]
            # insert new consensus
            val = (fixed_locus, coltag, rowtag, vseg, jseg, experiment_id)
            ins_consensus = "INSERT IGNORE INTO {}.consensus_stats (locus, col_tag, row_tag, best_V, best_J, experiment_id) VALUES (%s,%s,%s,%s,%s,%s) ON DUPLICATE KEY UPDATE consensus_id = LAST_INSERT_ID(consensus_id);".format(config_scireptor['database'])
            cursor.execute(ins_consensus,val)
            count_seq_ids+=1
            consensus_id = int(cursor.lastrowid)
                 # select the corresponding seq_ids
            sel_seqids = "SELECT reads.seq_id FROM %s.reads \
                          JOIN %s.reads_VDJ_segments as Vseg ON reads.seq_id = Vseg.seq_id \
                          JOIN %s.reads_VDJ_segments as Jseg ON reads.seq_id = Jseg.seq_id \
                          JOIN %s.sequencing_run ON reads.sequencing_run_id = sequencing_run.sequencing_run_id \
                          WHERE Vseg.type='V' AND Jseg.type='J'\
                          AND Vseg.igblast_rank=1 AND Jseg.igblast_rank=1 \
                          AND reads.well_id = '%s'\
                          AND sequencing_run.experiment_id = '%s' \
                          AND Vseg.VDJ_id = '%s' AND Jseg.VDJ_id = '%s' \
                          AND reads.consensus_id IS NULL;" %(config_scireptor['database'],config_scireptor['database'],config_scireptor['database'],config_scireptor['database'],well_id,experiment_id,vseg,jseg)
            cursor.execute(sel_seqids)
            n_seq=int(cursor.rowcount)
            sel_seqids_fetch = cursor.fetchall()
            update_nseqs = "UPDATE {}.consensus_stats SET n_seq = (n_seq+%d) WHERE consensus_id = %d;".format(config_scireptor['database'])
            cursor.execute(update_nseqs % (n_seq,consensus_id))
            sel_seqid =  []
            for seqid_row in sel_seqids_fetch:
                sel_seqid.append(seqid_row[0])
                cnt_reads_consensus_updated = 0
                if not sel_seqid:
                    pass
                else:
                    sql_ins_seq = "UPDATE %s.reads SET consensus_id=%s WHERE seq_id IN (%s);" % (config_scireptor['database'], consensus_id, ','.join(str(id) for id in sel_seqid))
                    cursor.execute(sql_ins_seq)
    print("Total number of well_ids processed:", count_wellids, "\n")
    print("Total number of assigned sequence ids:", count_seq_ids, "\n")
    print("[todb_consensus_tags.pl][INFO] Assigned", count_seq_ids, "sequence IDs to", count_wellids, "well IDs for locus", fixed_locus,"\n")
