#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 17/6/2024 12:39
# @Author  : Runsheng
# @File    : test_motif.py


#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2/13/2021 1:46 PM
# @Author  : Runsheng
# @File    : clustercj_test.py
"""
Test set for clustercj mode
"""
import os
from hammermotif.motif import *
import unittest


class HammermotifTest(unittest.TestCase):
    def setUp(self):

        wkdir="/data/Dropbox/temp/hammermotif_test"
        print(wkdir)
        print(os.listdir(wkdir))
        self.wkdir=wkdir
        self.bed_file_path = os.path.join(wkdir, 'potential_modification_site.bed')
        self.genome= os.path.join(wkdir, 'ecoli.fa')

    def test_motif(self):
        bed_file_path= self.bed_file_path
        genome_fasta_path = self.genome
        output_dir = os.path.join(self.wkdir, "motif_out")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        motif_enrichment_pipeline(
            self.bed_file_path,
            self.genome,
            motifs=None,  # Known motifs to include
            min_len=3,
            max_len=5,
            fdr_threshold=0.01
        )


    def test_ccwggg(self):
        pass







    def tearDown(self):
        self=None


