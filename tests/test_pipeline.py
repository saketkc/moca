"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

import os
import shutil
import unittest
from Bio import SeqIO
from moca.pipeline import Pipeline
from moca.bedoperations import fimo_to_sites

class TestPipeline(unittest.TestCase):
    """Test pipeline """

    def setUp(self):
        """Setup"""
        self.configuration_file = 'tests/data/application.cfg'
        self.fasta = 'tests/data/macsPeak.out.fasta'
        self.pipeline = Pipeline(self.configuration_file)

    def test_meme(self):
        """Test meme runner"""
        if os.path.exists('tests/data/out/meme_analysis'):
            shutil.rmtree('tests/data/out/meme_analysis')
        output = self.pipeline.run_meme(fasta_in=self.fasta, out_dir='tests/data/out/meme_analysis')
        #TODO Check if meme.txt is same and created
        assert output['exitcode'] == 0
        if os.path.exists('tests/data/out/meme_analysis'):
            shutil.rmtree('tests/data/out/meme_analysis')

    def test_fimo(self):
        """Test fimo runner"""
        fimo_output = self.pipeline.run_fimo(motif_file='tests/data/meme_analysis/meme.txt',
                                             sequence_file=self.fasta, out_dir=None, strargs=None)
        assert fimo_output['exitcode'] == 0

    def test_fimo_to_sites(self):
        fimo_file = 'tests/data/fimo_analysis/fimo.txt'
        record_dict = SeqIO.to_dict(SeqIO.parse(open(self.fasta), "fasta"))
        fimo_df = fimo_to_sites(fimo_file)
        for i, row in fimo_df.iterrows():
            record_id = row['sequence name']
            fimo_sequence = row['matched sequence']
            start = row['start']-1 # 0-based start
            end = row['stop'] # 1-based end
            strand = row['strand']
            fasta_sequence = record_dict[record_id].seq
            if strand == '+':
               assert str(fasta_sequence)[start:end] == fimo_sequence
            else:
               assert str(fasta_sequence[start:end].reverse_complement()) == fimo_sequence


