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
from hammermotif.motif_meme import *
import unittest


class HammermotifTest(unittest.TestCase):
    def setUp(self):

        wkdir="/home/li/myapp/Hammerhead-motif/test"
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

        results=main(
                self.bed_file_path,
                self.genome
            )
        print(results)


    def test_class_motif(self):
        from Bio import SeqIO
        regions = parse_bed_file(self.bed_file_path,  flank_size=7)
        genome = SeqIO.to_dict(SeqIO.parse(self.genome,format="fasta"))
        sequences = extract_sequences(regions, genome)

        from Bio import SeqIO
        import json
        from datetime import datetime
        from pathlib import Path

        def save_motifs(motifs, output_dir: str):
            """
            Save motifs to JSON and text files.
            """
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Save detailed JSON output
            json_output = {
                'timestamp': datetime.utcnow().isoformat(),
                'motifs': [motif.to_dict() for motif in motifs]
            }

            with open(output_dir / 'motifs.json', 'w') as f:
                json.dump(json_output, f, indent=2)

            # Save simple text output
            with open(output_dir / 'motifs.txt', 'w') as f:
                f.write(f"# Motif analysis results - {datetime.utcnow().isoformat()}\n")
                f.write("# Format: Consensus\tIUPAC\tNum_sites\tE-value\n\n")

                for motif in motifs:
                    f.write(f"{motif.consensus}\t{motif.iupac}\t{motif.nsites}\t{motif.evalue:.2e}\n")


        # Initialize motif enrichment with parallel processing
        me = MotifEnrichment(min_width=4, max_width=7, n_processes=20)

        # Find enriched motifs
        enriched_motifs = me.find_enriched_motifs(sequences, n_motifs=20, min_sites=100)

        # Save results
        save_motifs(enriched_motifs, "output_motifs")

        # Print results
        print("\nEnriched Motifs:")
        for i, motif in enumerate(enriched_motifs, 1):
            print(f"\nMotif {i}:")
            print(f"Consensus: {motif.consensus}")
            print(f"IUPAC: {motif.iupac}")
            print(f"Number of sites: {motif.nsites}")
            print(f"E-value: {motif.evalue:.2e}")


    def test_ccwggg(self):
        pass




    def tearDown(self):
        self=None


