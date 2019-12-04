#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de
    -adapted from Katharina Imkeller-
    
Complete the donor, sample, sort and event tables for high-throughput experiments. To each scenario of donor, sample, etc... you assign a numerical identifier and can then specify the corresponding wells in the matrix.

Necessary input files (extended):
	-m	<experiment_id>_metainfo.tsv file, tab delimited (use template 48_48_ or 240_256_metainfo.xls->worksheet2 and store as tsv with tabs). The first column is the numerical identifier that will be used in the plate layout. The first two rows of the TSV will be ignored, since they contain the headers in the current spreadsheet.
	-p	<experiment_id>_plate.tsv file (use template 48_48_ or 240_256_metainfo.xls->worksheet1 and store as tsv with tabs). When parsing, first row and column are ignored, they contain row and col numbers. The other 48*48 or 240*256 cells contain the identifier that already appeared in metainfo.tsv to specify which well contains what.
	-pb	<experiment_id>_platebarcodes.tsv file with the plate barcode corresponding to each plate number.
	-exp	experiment_id

1. Get information on matrix (48_48 e.g.) and plate layout (384 well plates, nrows, ncols e.g.)
2. Prepare database insertion, selection and update statements. Notably in the sequences table the event_id will be updated. Duplicates on unique keys will be ignored, not overwritten.
3. Open the input files.
4. From the plate_layout.tsv for each identifier, remember in a hash of arrays the corresponding wells (identified by row-col position). This allows you afterwards to find all wells and corresponding sequences with a certain sample, sort, donor... scenario. The event_id can then easily be updated in the sequences table. Only sequences that do not yet have an assigned event_id are updated (prevents problems when several matrices are stored in one database). From platebarcodes store platenr-barcode relations into hash.
5. Go through the metainfo.tsv and try to consecutively insert donor, sample, sort. On the event level, go through all corresponding wells, insert into event table and then update all correspnding sequences without event_id.

LOGGED INFORMATION

- donor, sample, sort information
- how many events belonged to the donor,sample,sort combi and how many sequences where found
- corrected tags in case applicable

"""

import argparse
from argparse import RawTextHelpFormatter
import bcelldb_init as binit
import glob, math


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-p', '--plate_input', required=True, help="<experiment_id>_metainfo.tsv")
parser.add_argument('-m', '--metainfo_input', required=True, help="<experiment_id>_plate.tsv file")
parser.add_argument('-pb', '--plate_barcodes', required=True, help="<experiment_id>_platebarcodes.tsv")
parser.add_argument('-exp', '--experiment_id', required=True, help="experiment_id")
args = parser.parse_args()


#parse config and get db connection
conf = binit.get_config()
db = binit.connect()
cursor = db.cursor()

config_dict_loci = {"Ig": ["H", "K", "L"],
                    "TCR": ["B", "A"]}

loci_current = config_dict_loci[conf["igblast_seqtype"]]

### 1. Get information on matrix and plate layout (how many wells to look for)

col_num, row_num = conf["matrix"].split("_")
n_col_per_plate = int(conf["ncols_per_plate"])
n_row_per_plate = int(conf["nrows_per_plate"])

### 2. Prepare DB insertion statement

# prepare statements for donor

ins_donor = ("INSERT INTO {}.donor ("
             "donor_identifier, background_treatment, project, strain, add_donor_info,"
             " species_id) VALUES (%s,%s,%s,%s,%s,%s"
             ") ON DUPLICATE KEY UPDATE donor_id = LAST_INSERT_ID(donor_id)".format(conf['database']))

# sample
ins_sample = ("INSERT INTO {}.sample (tissue, sampling_date, add_sample_info, donor_id"
              ") VALUES (%s,%s,%s,%s) ON DUPLICATE KEY UPDATE sample_id=LAST_INSERT_ID("
              "sample_id)".format(conf['database']))
# sort
ins_sort = ("INSERT INTO {}.sort (antigen, population, sorting_date, add_sort_info, sample_id"
            ") VALUES (%s,%s,%s,%s,%s) ON DUPLICATE KEY UPDATE sort_id=LAST_INSERT_ID("
            "sort_id)".format(conf['database']))

# event
ins_event = ("INSERT INTO {}.event (well, plate, row, col, sort_id, plate_layout_id, plate_barcode"
             ") VALUES (%s,%s,%s,%s,%s,%s,%s)  ON DUPLICATE KEY UPDATE event_id=LAST_INSERT_ID("
             "event_id)".format(conf['database']))

# select sequence_id where event_id will be inserted
sel_seq_id = ("SELECT seq_id FROM {}.sequences JOIN {"
              "}.consensus_stats ON consensus_stats.sequences_seq_id = sequences.seq_id \
              WHERE event_id IS NULL AND row_tag='%s' and col_tag='%s' AND sequences.locus='%s' AND experiment_id='%s'".format(
    conf['database'], conf['database']))

# update event_id
update_event = ("UPDATE {}.sequences SET event_id=%s where seq_id=%s".format(conf['database']))


### 3. Open infiles

plate = open(glob.glob(args.plate_input)[0], "r")
meta = open(glob.glob(args.metainfo_input)[0], "r")
barcodes = open(glob.glob(args.plate_barcodes)[0], "r")

id_r_c_dic = {}
count_line = 0

for line in plate:
    line = line.strip()
    count_line += 1
    if count_line > 1 or not count_line < int(row_num) + 1:
        row = count_line - 1
        entries = line.split("\t")
        # get entries from 2nd column on (first column contains index)
        for i in range(1,len(entries)):

            id_r_c_dic["R{}-C{}".format(str(row).zfill(3), str(i).zfill(3))]= entries[i]

# swap values as keys of dictionary, there should be a smarter way to do this already in that ^ loop
id_row_col_dic = {n:[k for k in id_r_c_dic.keys() if id_r_c_dic[k] == n] for n in set(id_r_c_dic.values())}

# print (id_row_col_dic)
plate.close()


### Extract plate barcode
plate_barcode_dic={}
count_line = 0

for line in barcodes:
    line=line.strip()
    count_line += 1
    if count_line > 1 and len(line.split("\t"))==2:
        plate_nr,plate_barcode = line.split("\t")
        plate_barcode_dic[plate_nr]=plate_barcode
        if plate_barcode and int(conf["log_level"]) >= 4:
            print("[todb_sampleinfo_highth.pl][DEBUG] Read plate barcode \"{}\" for plate {} from metadata file".format(plate_barcode, plate_nr))
        elif not plate_barcode and int(conf["log_level"]) >= 4:
            print ('[todb_sampleinfo_highth.py][DEBUG]no barcode provided for plate: {} '.format(plate_nr))


#### 5. Extract sample information for each id from METAqqq

count_line = 0

for line in meta:

    line = line.strip()
    count_line += 1

    if count_line > 2 and len (line) >= 10:

        (identifier, antigen,
        population, sorting_date, add_sort_info, tissue, sampling_date, add_sample_info,
        donor_identifier, background_treatment, project, strain, add_donor_info) = line.split("\t")
        
        # insert donor
        cursor.execute(ins_donor, (donor_identifier, background_treatment, project, strain, add_donor_info, conf["species"]))
        donor_id = cursor.lastrowid

        # insert sample
        cursor.execute(ins_sample,(tissue, sampling_date, add_sample_info, donor_id))
        sample_id = cursor.lastrowid

        # insert sort
        cursor.execute(ins_sort,(antigen, population, sorting_date, add_sort_info, sample_id))
        sort_id = cursor.lastrowid

        # count events and sequences for that donor-sample-sort combi
        count_events = 0
        count_sequences =0

        if id_row_col_dic[identifier]:
            wells = id_row_col_dic[identifier]
            if int(conf["log_level"]) >= 4:
                print("[todb_sampleinfo_highth.py][DEBUG] Metainfo identifier {} has {} wells on plates. ".format(identifier, str(len(wells))))

            for well in wells:
                row_tag, col_tag = well.split("-")
                row = int(row_tag[1:])
                col = int(col_tag[1:])

                # only if correct row and col have been found
                if int(row) > 0 and  int(col)>0:

                    # convert row col information to well plate
                    plate = math.ceil(col/int(
                        n_col_per_plate)) + ((math.ceil(row/n_row_per_plate - 1) * (
                            int(col_num)/n_col_per_plate)))

                    # for modulo calculation origin needs to be (0,0)
                    new_row = row - 1
                    new_col = col - 1
                    well = (new_col % n_col_per_plate) + (new_row % n_row_per_plate) * n_col_per_plate + 1

                    cursor.execute(ins_event, (well, plate, row, col,sort_id,
                                               conf['plate_layout'], plate_barcode_dic['{}'.format(int(plate))]))

                    count_events += 1
                    event_id = cursor.lastrowid

                    for locus in loci_current:

                        # todo correction rewrite # correct_tagconfusion.pl
                        corr_col_tag = col_tag
                        corr_row_tag = row_tag
                        # get the corresponding sequence id
                        cursor.execute(sel_seq_id % (corr_row_tag, corr_col_tag, locus, args.experiment_id))
                        row = cursor.fetchall() # added to fetch doublets of H/K/L on single event_id's

                        if row:
                            for each_seq in row:
                                seq_id = each_seq[0]
                                cursor.execute(update_event % (event_id, seq_id))
                                count_sequences += 1

            if int(conf["log_level"]) >= 3:

                print("[todb_sampleinfo_highth.py][INFO] Donor \"{}\", sample \"{}\", sort \"{}\": {} events and {} sequences ( Loci: {} ).".format(
                    donor_id, sample_id, sort_id, count_events, count_sequences, loci_current))

        else:

            print("[todb_sampleinfo_highth.py][WARNING] Metadata identifier {} (Donor ID: {} Sample ID: {} Sort ID: {}) has NO assigned events in _plate.tsv.".format(identifier, donor_id, sample_id, sort_id))