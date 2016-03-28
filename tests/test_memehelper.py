"""
test_motifopearations
----------------------------------

Tests for `moca.helpers.meme` module.
"""

import os
import shutil
import unittest
from Bio import SeqIO
from moca.helpers import read_memefile
from moca.helpers import get_motif_ic
from moca.helpers import get_motif_bg_freq
from moca.helpers import get_total_sequences
import numpy as np

class TestMemeOperations(unittest.TestCase):
    """Test pipeline """

    def setUp(self):
        """Setup"""
        self.meme_file = 'tests/data/expected_out/meme_analysis/meme.txt'

    def test_memeprofile(self):
        """ Test entropy caclulations
	         bits    2.3
                 2.1
                 1.8                        *
                 1.6                        **
Relative         1.4                *     * **
Entropy          1.2 *      **  *** *     * **
(23.5 bits)      0.9 ** * * **  *******  ** **
                 0.7 **** * **  *******  *****
                 0.5 ********** **************
                 0.2 *************************
                 0.0 -------------------------


        """
        motifs = read_memefile(self.meme_file)
        record = motifs['motif_records'][0]
        motif_ic = get_motif_ic(self.meme_file, 0)
        target = np.array([1.2,0.9,0.7,0.9,0.5,
                           0.9,0.5,1.2,1.2,0.5,
                           0.2,1.2,1.2,1.2,0.9,
                           1.4,1.2,1.2,0.5,0.9,
                           1.2,1.4,0.7,1.8,1.6])

        assert np.allclose(target, motif_ic, atol=0.6)

    def test_bgfreq(self):
        bg_freq = get_motif_bg_freq(self.meme_file)
        assert bg_freq == {'A':0.202, 'C':0.298, 'G': 0.298, 'T': 0.202}

    def test_totalseq(self):
        assert get_total_sequences(self.meme_file) == 10
