import os
import re
from moca.helpers import MocaException
import numpy as np

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
        third_line = f.readline().strip()
        second_line = re.sub( '\s+\t+', ' ', second_line)
        third_line = re.sub( '\s+\t+', ' ', third_line)

        columns = second_line.split()
        values = third_line.split()

        if len(columns) != len(values):
            raise MocaException('Error parsing centrimo file. \
                                This is a bug and should be reported upstream')

        return dict(zip(columns, values))


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
    data = np.genfromtxt(os.path.abspath(centrimo_stats),
                         unpack=True,
                         skip_header=True)
    return data[0,:], data[1,:]

