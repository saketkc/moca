"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

import unittest
import filecmp
from moca.bedoperations import Bedfile
from moca.helpers import MocaException

class TestBedoperations(unittest.TestCase):

    def setUp(self):
        self.narrowpeak = 'tests/data/narrowPeak.bed'
        self.broadpeak = 'tests/data/broadPeak.bed'
        self.macspeak = 'tests/data/macsPeak.bed'
        self.getfastapeak = 'tests/data/getfasta.bed'
        self.genome_table = 'tests/data/hg19.chrom.sizes'

    def test_narrowPeak(self):
        """Test load narrowPeak"""
        loaded_bed = Bedfile(self.narrowpeak, self.genome_table)
        assert loaded_bed.bed_format == 'narrowPeak'

    def test_broadPeak(self):
        """Test load broadPeak"""
        with self.assertRaises(MocaException):
            Bedfile(self.broadpeak, self.genome_table)

    def test_macsPeak(self):
        """Test load macsPeak"""
        loaded_bed = Bedfile(self.macspeak, self.genome_table)
        assert loaded_bed.bed_format == 'macsPeak'

    def test_generatefasta(self):
        """Test generate fasta"""
        loaded_bed = Bedfile(self.getfastapeak, 'tests/data/getfasta.gt')
        #loaded_bed.determine_peaks()
        t = loaded_bed.slop_bed(flank_length=20)
        print t#loaded_bed
        fasta_out = loaded_bed.extract_fasta(fasta_in='tests/data/getfasta.fa', fasta_out='tests/data/generated_out/getfasta.fa')
        assert filecmp.cmp('tests/data/generated_out/getfasta.fa', 'tests/data/expected_out/getfasta.expected.fa')

    def test_scorefile(self):
        loaded_bed = Bedfile(self.narrowpeak, self.genome_table)
        assert filecmp.cmp('tests/data/narrowPeak.sorted','tests/data/expected_out/narrowPeak.sorted')
