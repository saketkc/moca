"""
test_wigoperations
----------------------------------

Tests for `moca.wigoperations` module.
"""

import unittest
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
        assert scores == [[0.10000000149011612, 0.20000000298023224, 0.30000001192092896].reverse(),
                          [1.5]]


