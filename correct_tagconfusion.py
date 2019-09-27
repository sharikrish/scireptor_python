#!/usr/bin/python3

"""
@franasa - f.arcila@dkfz-heidelberg.de
    -adapted from Katharina Imkeller-

correct_tagconfusion

This module is called by todb_sampleinfo_highth.py in order to correctly assign the event_id, even when when tags are in the wrong positions (problems
in the primer matrix) or plates got mixed up during processing.

Use correct_tags($tag_column, $tag_row, $locus)) with $tag_(column|row) =~ [RC][0-9][0-9][0-9] and $locus =~ [HKL] to return correct tag name. The
information which of the following correction to apply is set via the "tag_correction" key in the config file:

B< none >	(implemented Jan 2014): No correction (default).

B< tag_batch >	(implemented Jan 2014):	Correct for a switch of lambda tags in the tag batch Mm_AGM_002 240x256 matrix. This switch occurred during
primer manufacturing and is already	present	in the primer master plate:
For plate rows 1-12 (physical rows 1-192), lambda reads located in the even physical rows and odd plate rows appear to be located one plate row
(16 physical rows) below, while those located in even physical rows and even plate rows appear to be located one plate row above. Reads in odd physical
rows appear at the correct position.

B< D01_plateposition >	(implemented Feb 2014):	Correct for plate swapping in D01 experiment:
Kappa chain Row1-16:Col25-48 is exchanged by Row17-32:Col1-24 and vice versa.

"""


#### mouse matrix 240_256, tags mixed due partial quadrant exchange during synthesis

class CorrectConfusions:

    def __init__(self,row_tag, col_tag, locus):
        self.breed = breed
        self.row_tag = row_tag
        self.col_tag = col_tag
        self.locus = locus


    def correct_tags_240_256_1 (self, old_tag_1):

        # Kappa chain Row1-16:Col25-48 is exchanged by Row17-32:Col1:24 and vice versa

        self.old_col_tag = old_tag_1
        old_row_tag = ""




    ### Experiment D01

    def correct_tags_D01_2():
        pass
