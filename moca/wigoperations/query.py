"""Perform query on bigwig files
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object
import warnings
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
        chrom_lengths = self.get_chromosomes
        for chrom, chromStart, chromEnd, strand in intervals:
            if chrom not in list(chrom_lengths.keys()):
                warnings.warn('Chromosome {} does not appear in the bigwig'.format(chrom), UserWarning)
                continue

            chrom_length = chrom_lengths[chrom]
            if int(chromStart)> chrom_length:
                raise MocaException('Chromsome start point exceeds chromosome length: {}>{}'.format(chromStart, chrom_length))
            elif int(chromEnd)> chrom_length:
                raise MocaException('Chromsome end point exceeds chromosome length: {}>{}'.format(chromEnd, chrom_length))
            score = self.wig.values(chrom, int(chromStart), int(chromEnd))
            if strand == '-':
                score.reverse()
            scores.append(score)
        return np.array(scores)

    @property
    def get_chromosomes(self):
        """Return list of chromsome and their sizes
        as in the wig file

        Returns
        -------
        chroms: dict
            Dictionary with {"chr": "Length"} format

        """
        return self.wig.chroms()
