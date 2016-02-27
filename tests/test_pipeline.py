"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

import os
import shutil
import unittest
from moca.pipeline import Pipeline

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
        self.pipeline.run_meme(fasta_in=self.fasta, out_dir='tests/data/out/meme_analysis')
        fimo_output = self.pipeline.run_fimo(motif_file='tests/data/out/meme_analysis/meme.txt',
                                             sequence_file=self.fasta, out_dir=None, strargs=None)
        assert fimo_output['exitcode'] == 0

