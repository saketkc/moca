"""
test_wigoperations
----------------------------------

Tests for `moca.wigoperations` module.
"""

import unittest
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


