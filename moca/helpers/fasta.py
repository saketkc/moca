import os
from Bio import motifs
import glob
import logging
import sys
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import seed
from random import randint

def make_uppercase(mixed_fasta, upper_fasta):
    """Convert fasta to have all upper case letters

    This function is essentially a hack to convert
    repeat masked sequences which are often repeated by MEME
    to normal alphabets. The problem arises
    becauses Biopython handles only upper case letters(IUPACUnambiguousDNA) in the
    core motif model. In cases where MEME reports motifs from
    soft-masked regions, the meme file can't be read by Biopython
    unless it is generated from a fasta containing only upper case
    letters. MEME's recommendation is to NOT use repeat masked
    sequence. It is a bit surprising in the first place thatMEME reports regions
    from repeatmasked sites.
    """
    records = (rec.upper() for rec in SeqIO.parse(os.path.abspath(mixed_fasta), 'fasta'))
    SeqIO.write(records, os.path.abspath(upper_fasta), 'fasta')

def generate_random_fasta(genome, genome_table, random_fasta):
    """Generate fasta pooling seqe"""
    gt = None
    seed(1)
    with open(genome_table) as f:
        gt = f.readlines()
    chr_map = {}
    for line in gt:
        line = line.strip()
        s = line.split('\t')
        chr_map[s[0]] = int(s[1])
    records = list(SeqIO.parse(open(genome,'r'), 'fasta'))
    seqs = []
    chr_keys = chr_map.keys()
    for i in range(0,MAX_PEAKS_TO_KEEP):
        ## Avoid regions from scaffolds
        chr_idx = None
        chr_id = None
        while (chr_idx is None):
            chr_index = randint(0, len(chr_keys)-1)
            chr_id = chr_keys[chr_index]
            split = chr_id.split('_')
            if len(split)==1:
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
        start = randint(0, limit-FLANKING_SEQ_LENGTH-1)
        end = start+FLANKING_SEQ_LENGTH
        data = record.seq[start:end]
        seq = SeqRecord(data,'{}_{}_-{}'.format(chr_id, start, MOTIF_FLANKING_BASES),'','')
        seqs.append(seq)

    output_handle = open(random_fasta, 'w')
    logging.info('###########GeneratingRandomFA Start########################')
    SeqIO.write(seqs, output_handle, 'fasta')
    logging.info('###########GeneratingRandomFA End########################')
    output_handle.close()
    return random_fasta

