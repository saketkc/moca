"""Fasta helper functions"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import os
from Bio import SeqIO
import pandas
from pyfaidx import Fasta
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
    length = len(record_dict[record_dict.keys()[0]])
    for key in record_dict.keys():
        assert length == len(record_dict[key])
    fasta_metadata = {'num_seq': len(record_dict.keys()), 'len_seq': length}
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

def generate_random_fasta(genome,
                       genome_table,
                       num_seq,
                       len_seq,
                       out_fasta):
    """Generate fasta pooling seqe
    Attributes
    ---------
    genome: str
        Path to genome
    genome_table: str
        Path to chromosome size
    len_seq: int
        Length of fasta sequences to generate
    num_seq: int
        Number of fasta records to generate
    out_fasta: str
        Path to write random fasta
    """
    filt_func = lambda chrom: '_' not in chrom and chrom[-1].isdigit()
    gt_map = pandas.read_table(genome_table, index_col=0, header=None)
    chr_keys = gt_map.index.tolist()
    ## Avoid scaffolds and enforce the last chracter of chromsome key being numeric
    ## TODO This is not  foolproofi. though I can't think of cases it will fail, but it is also not random in true sense
    chr_keys_filtered = filter(lambda chrom: '_' not in chrom and chrom[-1].isdigit(), chr_keys)
    chr_selected = np.random.choice(chr_keys_filtered, num_seq)

    chr_selected_length = gt_map.ix[chr_selected].values.flatten()
    chr_selected_start = np.array([np.random.choice(np.arange(1, chr_length-len_seq)) for chr_length in chr_selected_length])
    chr_selected_end = chr_selected_start+len_seq-1

    fasta = Fasta(genome, read_ahead=10000, filt_function=filt_func, sequence_always_upper=True)

    with open(out_fasta, 'w') as f:
        for chr_name, chr_start, chr_end in zip(chr_selected, chr_selected_start, chr_selected_end):
            seq = fasta[chr_name][chr_start:chr_end]
            f.write('>{}:{}-{}\n'.format(seq.name, seq.start, seq.end))
            f.write(seq.seq + '\n')
