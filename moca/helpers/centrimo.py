from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import zip
import os
import re
from ..helpers import MocaException

def read_centrimo_txt(centrimo_txt):
    """Read centrimo.txt for single motif

    Attributes
    ----------
    centrimo_txt: str
        Path to centrimo.txt file

    Returns
    -------
    hits: dict
         A dict containing hits from centrimo.txt

    """
    with open(os.path.abspath(centrimo_txt)) as f:
        assert f.readline().strip() == '# WARNING: this file is not sorted!'
        second_line = f.readline().strip().replace('# ','')
        second_line = re.sub( '\s+\t+', ' ', second_line)
        values = []
        should_read = True
        while should_read:
            value = f.readline().strip()
            if value == '##':
                break
            value = re.sub('\s+\t+', ' ', value)
            values.append(value)
        columns = second_line.split()
        if len(columns) != len(values[0].split()):
            raise MocaException('Error parsing centrimo file. \
                                This is a bug and should be reported upstream')
        centrimo_dicts = [dict(list(zip(columns, val.split()))) for val in values]
        return centrimo_dicts


def read_centrimo_stats(centrimo_stats):
    """Read centrimo' stats.txt for single motif[should be independent though]

    Attributes
    ----------
    centrimo_stats: str
        Path to stats.txt file

    Returns
    -------
    hits: dict
         A dict containing hits from centrimo.txt

    """
    motif_wise_site_stats = {'MEME': {}, 'DREME': {}}
    with open(centrimo_stats) as f:
        for line in f:
            line = line.strip()
            if line.startswith('DB'):
                line_split = line.replace(' ', '\t').split('\t')
                assert len(line_split) == 5
                meme_dreme = line_split[-1]
                assert meme_dreme in ['MEME', 'DREME']
                key = ('_').join(line_split[2:4])
                assert key not in list(motif_wise_site_stats[meme_dreme].keys())
                motif_wise_site_stats[meme_dreme][key] = {'pos': [], 'count': []}
            else:
                line = line.replace(' ', '')
                pos, count = line.split('\t')
                motif_wise_site_stats[meme_dreme][key]['pos'].append(float(pos))
                motif_wise_site_stats[meme_dreme][key]['count'].append(float(count))

    return motif_wise_site_stats

