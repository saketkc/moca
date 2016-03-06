"""Fasta helper functions"""
import logging
import os
from random import seed
from random import randint
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
    #TODO Make this a ppandas indexed dataframe?
    gt_map = None
    seed(1)
    with open(os.path.abspath(genome_table)) as f:
        gt_map = f.readlines()
    chr_map = {}
    for line in gt_map:
        line = line.strip()
        s = line.split('\t')
        chr_map[s[0]] = int(s[1])
    records = list(SeqIO.parse(open(os.path.abspath(genome), 'r'), 'fasta'))
    seqs = []
    chr_keys = chr_map.keys()
    for i in range(0, num_seq):
        ## Avoid regions from scaffolds
        chr_idx = None
        chr_id = None
        while chr_idx is None:
            chr_index = randint(0, len(chr_keys)-1)
            chr_id = chr_keys[chr_index]
            split = chr_id.split('_')
            if len(split) == 1:
                chr_idx = chr_id
            else:
                pass
        record = None
        if chr_id is None:
            raise RuntimeError('Chr_id cannot be none')

        for r in records:
            if r.id == chr_id: #'chr{}'.format(chr_index):
                record = r#records[chr_index]
        if not record:
            raise RuntimeError('No matching chrosome found for: {}'.format(chr_id))

        ##We use limit from genome_table!
        limit = chr_map[chr_id]
        start = randint(0, limit-len_seq-1)
        end = start+len_seq-1
        data = record.seq[start:end]
        seq = SeqRecord(data, '{}:{}-{}'.format(chr_id, start, end), '', '')
        seqs.append(seq.upper())

    output_handle = open(out_fasta, 'w')
    logging.info('###########GeneratingRandomFA Start########################')
    SeqIO.write(seqs, output_handle, 'fasta')
    logging.info('###########GeneratingRandomFA End########################')
    output_handle.close()

