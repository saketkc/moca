"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

from moca.bedoperations import Bedfile
import unittest
from moca.helpers import MocaException

class TestBedoperations(unittest.TestCase):

    def setUp(self):
        self.narrowpeak = 'tests/data/narrowPeak.bed'
        self.broadpeak = 'tests/data/broadPeak.bed'
        self.macspeak = 'tests/data/macsPeak.bed'
        self.genome_table = 'tests/data/hg19.chrom.sizes'

    def test_narrowPeak(self):
        loaded_bed = Bedfile(self.narrowpeak, self.genome_table)
        assert loaded_bed.bed_format == 'narrowPeak'

    def test_broadPeak(self):
        with self.assertRaises(MocaException):
            loaded_bed = Bedfile(self.broadpeak, self.genome_table)

    def test_macsPeak(self):
        loaded_bed = Bedfile(self.macspeak, self.genome_table)
        assert loaded_bed.bed_format == 'macsPeak'

