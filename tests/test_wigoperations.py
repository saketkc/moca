"""
test_wigoperations
----------------------------------

Tests for `moca.wigoperations` module.
"""

from moca.wigoperations import WigReader
import unittest
from moca.helpers import MocaException
from array import array
from math import isnan

class TestWigoperations(unittest.TestCase):

    def setUp(self):
        self.wig_location = 'tests/data/test.bw'

    def test_open(self):
        loaded_wig = WigReader(self.wig_location)
        assert(loaded_wig is not None)

    def test_query(self):
        loaded_wig = WigReader(self.wig_location)
        scores = loaded_wig.query([("1", 0, 3, '-'),("1",150,151,'+')])
        assert scores == [[0.10000000149011612, 0.20000000298023224, 0.30000001192092896].reverse(), [1.5]]


