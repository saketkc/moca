from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import map
from builtins import zip
from builtins import range
import os
import operator
from math import log
from Bio import motifs
import numpy as np

def get_motif_bg_freq(meme_file):
    """Read backgroud frequencies as determined by MEME
    Biopython does not support reading background frequencies from MEME file
    This method implements the missing feature.

    Parameters
    ----------
    meme_file: str
        Location of MEME file

    Returns
    -------
    bg_frequencies = dict
        A dict with bases as keys and corresponding frequency  as values
    """
    frequencies_line = None
    with open(meme_file) as f:
        for line in f:
            if 'Background letter frequencies (from' in line.strip():
                frequencies_line = f.next().strip()

    assert frequencies_line is not None
    frequencies = frequencies_line.split(' ')
    values = [float(x) for x in frequencies[1:][::2]]
    keys = frequencies[0:][::2]
    bg_frequencies = {k:v for k,v in zip(keys, values)}
    return bg_frequencies

def get_total_sequences(meme_file):
    """Get total sequences
    Parameters
    ----------
    meme_file: str
        Location of MEME file

    Returns
    -------
    num_seq: int
        Number of sequences that meme was used on
    """
    seq_line = None
    count = 0
    with open(meme_file) as f:
        for line in f:
            if line.strip().startswith('TRAINING SET'):
                star_line = f.next().strip()
                datafile_line = f.next().strip()
                alphabet_line = f.next().strip()
                header_line = f.next().strip()
                dashed_line = f.next().strip()
                count = 0
                while not f.next().strip().startswith('*********'):
                    count+=1
    return count*2

#TODO Rename this!
def read_memefile(meme_file):
    """Summariser for MEME file
    Read meme file
    Parameters
    ----------
    meme_file: str
        Location of MEME file

    Returns
    -------
    summary: dict
        A summary containing the following details:
            - motif_occurences: dict contatining number of type each motif occurs. dict is indexed by key: 'motif1', 'motif2' etc
            - motif_records: List of Biopython motif objects
    summary:
    """
    summary = {}
    summary['motif_occurrences'] = {}
    records = motifs.parse(open(os.path.abspath(meme_file)), 'meme')
    summary['total_motifs'] = len(records)
    num_occurrences = []
    for index, record in enumerate(records):
        num_occurrences.append(int(getattr(record,'num_occurrences','Unknown')))

    sorted_occurences = sorted(enumerate(num_occurrences), key=lambda x: x[1])
    summary['motif_occurrences'] = {'motif{}'.format(index+1):value for index,value in sorted_occurences}
    summary['motif_records'] = records
    ### Read background frequenceies H since bioppython does not support them
    bg_frequencies = get_motif_bg_freq(meme_file)
    summary['bg_frequencies'] = bg_frequencies
    return summary

def create_position_profile(record, count_type='pwm'):
    """Find bases with minimum or maximum count type"""
    bases = ['A', 'C','T','G']
    base_profile = getattr(record, count_type)
    base_A = np.array(base_profile['A'])
    base_C = np.array(base_profile['C'])
    base_T = np.array(base_profile['T'])
    base_G = np.array(base_profile['G'])
    profile = np.vstack((base_A,base_C,base_T,base_G))
    position_profile = []
    for col in range(0, record.length):
        p_profile = profile[:,col]
        position_profile.append({k:v for k,v in zip(bases,p_profile)})
    return position_profile

def get_max_occuring_bases(record, max_count, count_type='counts'):
    """Find bases with maximum frequency at each position of motif record

    Given a motif record, find at each position, the base with maximum
    frequency and it's frequency

    """
    profile = create_position_profile(record, count_type)
    sorted_profile = []
    for i, p in enumerate(profile):
        sorted_p = sorted(list(p.items()), key=operator.itemgetter(1))
        sorted_profile.append(sorted_p[-max_count:])
    return sorted_profile

def position_wise_profile(counts_dict, length):
    """
    Convert base to position wise profile
    """
    profile = list(map(dict, list(zip(*[[(k, v) for v in value] for k, value in list(counts_dict.items())]))))
    return profile

def find_max_occurence(profile, max_count=2):
    """
    Return profile with base corresponding to max scores[CHECK!]
    """
    sorted_profile = []
    for p in profile:
        sorted_profile.append(sorted(list(p.items()), key=lambda x:x[1]))
    for i,p in enumerate(sorted_profile):
        sorted_profile[i] = p[-max_count:]
    return sorted_profile

def get_motif_ic(meme_file, n_motif):
    """Get position wise infomartion content of motif

    Parameters
    ----------
    motif: motif
        A Biopython motif record

    Returns
    -------
    motif_ic: array
        An array of information content at each base position
    """
    bases = ['A', 'C', 'T', 'G']
    motif_ic = []
    meme_summary = read_memefile(meme_file)
    motif = meme_summary['motif_records'][n_motif]
    bg_frequencies = meme_summary['bg_frequencies']
    log_odds = motif.pwm.log_odds()#bg_frequencies)
    pwm = motif.pwm
    for i in range(0, motif.length):
        s = 0
        for base in bases:
            p = pwm[base][i]
            if p<=0:
                s1 = 0
            else:
                s1 = log(p,2)
            s2 = log(bg_frequencies[base],2)
            s = s + p*(s1-s2)
        motif_ic.append(s)
    return motif_ic
