"""
test_wigoperations
----------------------------------

Tests for `moca.wigoperations` module.
"""

import unittest
import warnings
import numpy as np
from moca.wigoperations import WigReader

class TestWigoperations(unittest.TestCase):
    """Test Wigoperations"""

    def setUp(self):
        self.wig_location = 'tests/data/test.bw'

    def test_open(self):
        """Test load wig"""
        loaded_wig = WigReader(self.wig_location)
        assert loaded_wig is not None

    def test_query(self):
        """Test wig query"""
        loaded_wig = WigReader(self.wig_location)
        scores = loaded_wig.query([("1", 0, 3, '-'), ("1", 150, 151, '+')])
        expected_scores = np.array([[0.30000001192092896, 0.20000000298023224, 0.10000000149011612],
                          [1.5]])
        np.equal(expected_scores, scores)

    def test_chroms(self):
        loaded_wig = WigReader(self.wig_location)
        assert loaded_wig.get_chromosomes == {'1': 195471971, '10': 130694993}

    def test_out_of_bounds(self):
        """Test if query is out of bounds"""
        loaded_wig = WigReader(self.wig_location)
        with warnings.catch_warnings(record=True) as w:
            scores = loaded_wig.query([("M", 0, 3, '-')])
            assert len(w) >= 1
            assert str(w[-1].message) == 'Chromosome M does not appear in the bigwig'



