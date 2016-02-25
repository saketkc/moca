"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

from moca.pipeline import Pipeline
import os
import shutil
import unittest
from moca.helpers import MocaException

class TestPipeline(unittest.TestCase):

    def setUp(self):
        self.configuration_file = 'tests/data/application.cfg'
        self.fasta = 'tests/data/macsPeak.out.fasta'
        self.pipeline = Pipeline(self.configuration_file)

    def test_meme(self):
        if os.path.exists('tests/data/out/meme_analysis'):
            shutil.rmtree('tests/data/out/meme_analysis')
        output = self.pipeline.run_meme(fasta_in=self.fasta, out_dir='tests/data/out/meme_analysis')
        #TODO Check if meme.txt is same and created?
        assert output['exitcode'] != 1
        if os.path.exists('tests/data/out/meme_analysis'):
            shutil.rmtree('tests/data/out/meme_analysis')

