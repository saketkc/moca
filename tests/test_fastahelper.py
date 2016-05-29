"""
Tests for fasta helper
"""
import os
import shutil
import unittest
from moca.helpers import fasta

class TestFastaHelper(unittest.TestCase):
    """Test fasta helper"""
    def setUp(self):
        self.fasta = 'tests/data/getfasta.fa'

    def test_metadata(self):
        metadata = fasta.get_fasta_metadata(self.fasta)
        assert metadata['num_seq'] == 1
        assert metadata['len_seq'] == 205

