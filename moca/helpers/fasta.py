"""Fasta helper functions"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import zip
import os
from Bio import SeqIO
import pandas
import numpy as np

def get_fasta_metadata(in_fasta):
    """Return metadata about fasta
    Attributes
    ----------
    in_fasta: str
        Path to fasta file

    Returns
    -------
    fasta_metadata: dict
        A dict of containing the following fields: {'len_seq': length, 'num_seq': num}
    """
    record_dict = SeqIO.to_dict(SeqIO.parse(open(in_fasta), 'fasta'))
    length = len(record_dict[list(record_dict.keys())[0]])
    for key in list(record_dict.keys()):
        assert length == len(record_dict[key])
    fasta_metadata = {'num_seq': len(list(record_dict.keys())), 'len_seq': length}
    return fasta_metadata

def make_uppercase_fasta(mixed_fasta, upper_fasta):
    """Convert fasta to have all upper case letters

    This function is essentially a hack to convert
    repeat masked sequences which are often repeated by MEME
    to normal alphabets. The problem arises
    becauses Biopython handles only upper case letters(IUPACUnambiguousDNA) in the
    core motif model. In cases where MEME reports motifs from
    soft-masked regions, the meme file can't be read by Biopython
    unless it is generated from a fasta containing only upper case
    letters. MEME's recommendation is to NOT use repeat masked
    sequence.

    See: https://groups.google.com/forum/#!topic/meme-suite/0XEgRn0Lmcc

    Arguments
    ---------
    mixed_fasta: str
        Location of fasta with mixed(both uppercase and lowercase alphabets)
    upper_fasta: str
        Location to store the uppercase fasta
    """
    records = (rec.upper() for rec in SeqIO.parse(os.path.abspath(mixed_fasta), 'fasta'))
    SeqIO.write(records, os.path.abspath(upper_fasta), 'fasta')

