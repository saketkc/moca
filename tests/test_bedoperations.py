"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

import unittest
import filecmp
import pandas as pd
import numpy as np
from moca.bedoperations import Bedfile
from moca.helpers import MocaException

class TestBedoperations(unittest.TestCase):

    def setUp(self):
        self.narrowpeak = 'tests/data/narrowPeak.bed'
        self.broadpeak = 'tests/data/broadPeak.bed'
        self.macspeak = 'tests/data/macsPeak.bed'
        self.getfastapeak = 'tests/data/getfasta.bed'
        self.unsortedbed = 'tests/data/unsorted.bed'
        self.genome_table = 'tests/data/hg19.chrom.sizes'

    def test_narrowPeak(self):
        """Test load narrowPeak"""
        loaded_bed = Bedfile(self.narrowpeak, self.genome_table)
        assert loaded_bed.bed_format == 'narrowPeak'

    def test_broadPeak(self):
        """Test load broadPeak"""
        with self.assertRaises(MocaException):
            loaded_bed = Bedfile(self.broadpeak, self.genome_table)
            assert loaded_bed.bed_format == 'broadPeak'

    def test_macsPeak(self):
        """Test load macsPeak"""
        loaded_bed = Bedfile(self.macspeak, self.genome_table)
        assert loaded_bed.bed_format == 'macsPeak'

    def test_generatefasta(self):
        """Test generate fasta"""
        loaded_bed = Bedfile(self.getfastapeak, 'tests/data/getfasta.gt')
        total_peaks = loaded_bed.get_total_peaks
        train_bed_file, test_bed_file = loaded_bed.split_train_test_bed(train_peaks_count=total_peaks/2.0,
                                                                        test_peaks_count=total_peaks/2.0)

        train_slopped = loaded_bed.slop_bed(train_bed_file, flank_length=20)
        #test_sloppped = loaded_bed.slop_bed(test_bed_file, flank_length=20)

        train_fasta_out = loaded_bed.extract_fasta(bed_in=train_slopped,
                                                   fasta_in='tests/data/getfasta.fa',
                                                   fasta_out='tests/data/generated_out/getfasta.fa')
        assert filecmp.cmp('tests/data/generated_out/getfasta.fa', 'tests/data/expected_out/getfasta.fa')

    def test_slopper(self):
        loaded_bed = Bedfile(self.unsortedbed, 'tests/data/getfasta.gt')
        total_peaks = loaded_bed.get_total_peaks
        train_bed_file, test_bed_file = loaded_bed.split_train_test_bed(train_peaks_count=total_peaks/2.0,
                                                                        test_peaks_count=total_peaks/2.0)
        assert filecmp.cmp('tests/data/unsorted.train', 'tests/data/expected_out/unsorted.train')

    def test_scorefile(self):
        loaded_bed = Bedfile(self.narrowpeak, self.genome_table, 'tests/data/generated_out')
        assert filecmp.cmp('tests/data/narrowPeak.sorted','tests/data/expected_out/narrowPeak.sorted')

    def test_sort(self):
        bedfile = 'tests/data/ENCSR000AKB/ENCFF002CDP_1000_unsorted.bed'
        loaded_bed = Bedfile(bedfile, self.genome_table, 'tests/data/generated_out')
        df1 = pd.read_table('tests/data/generated_out/ENCFF002CDP_1000_unsorted.sorted')
        df2 = pd.read_table('tests/data/ENCSR000AKB/ENCFF002CDP_1000.bed')
        assert np.all(df1[[0]]==df2[[0]])
        assert np.all(df1[[4]]==df2[[4]])


    def test_traintest(self):
        bedfile = 'tests/data/ENCSR000AKB/ENCFF002CDP_1000_unsorted.bed'
        loaded_bed = Bedfile(bedfile, self.genome_table, 'tests/data/generated_out')
        total_peaks = loaded_bed.get_total_peaks
        train_bed_file, test_bed_file = loaded_bed.split_train_test_bed(train_peaks_count=total_peaks/2.0,
                                                                        test_peaks_count=total_peaks/2.0)
        assert filecmp.cmp('tests/data/generated_out/ENCFF002CDP_1000_unsorted.train',
                           'tests/data/expected_out/ENCFF002CDP_1000_unsorted.train')

