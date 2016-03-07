"""Perform query on bigwig files
"""
import numpy as np
import pyBigWig
from ..helpers import MocaException

class WigReader(object):
    """Class for reading and querying wigfiles"""
    def __init__(self, wig_location):
        """
        Arguments
        ---------
        wig_location: Path to bigwig

        """
        self.wig_location = wig_location
        try:
            self.wig = pyBigWig.open(self.wig_location)
        except Exception as e:
            raise MocaException('Error reading wig file: {}'.format(e))

    def query(self, intervals):
        """ Query regions for scores

        Arguments
        ---------
        intervals: list of tuples
            A list of tuples with the following format: (chr, chrStart, chrEnd, strand)

        Returns
        -------
        scores: np.array
            A numpy array containing scores for each tuple
        """
        scores = []
        for chrom, chromStart, chromEnd, strand in intervals:
            score = self.wig.values(chrom, int(chromStart), int(chromEnd))
            if strand == '-':
                score.reverse()
            scores.append(score)
        return np.array(scores)

    def _get_chromosomes(self):
        """Return list of chromsome and their sizes
        as in the wig file

        Returns
        -------
        chroms: dict
            Dictionary with {"chr": "Length"} format

        """
        return self.wig.chroms()
