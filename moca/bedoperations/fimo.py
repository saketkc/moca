"""Convert fimo.txt to a bed file"""

import os
import pandas as pd
from moca.helpers import filename_extension


def fimo_to_sites(fimo_file):
    """Convert fimo.txt to bed file
    fimi.txt columns:
    #pattern name:  The motif identifier
    #sequence name: The sequence identiifer
    strand:         The strand + indicates the motif matched the forward strand,\
                    - the reverse strand, and . indicates strand is\
                    not applicable (as for amino acid sequences).
    start:          The start position of the motif occurence (closed, 1-based coordinates,\
                    unless genomic coordinates are provided)\
    stop:           The end position of the motif occurence\
                    (closed, 1-based coordinates, unless genomic coordinates are provided).
    score:          The score for the motif occurence.\
                    The score is computed by by summing the appropriate entries from\
                    each column of the position-dependent scoring matrix that represents the motif.
    p-value:        The p-value of the motif occurence.\
                    The p-value is the probability of a random sequence of\
                    the same length as the motif matching that position of\
                    the sequence with a score at least as good.
    q-value:        The q-value of the motif occurence.\
                    The q-value is the estimated false discovery rate if the occurrence is\
                    accepted as significant.
                    See Storey JD, Tibshirani R. Statistical significance for genome-wide studies.\
                            Proc. Natl Acad. Sci. USA (2003) 100:9430-9445
    sequence:        The sequence matched to the motif.
    """
    fimo_file = os.path.abspath(fimo_file)
    fimo_df = pd.read_table(fimo_file)
    filename, ext = filename_extension(os.path.abspath(fimo_file))
    fimo_sites = os.path.join(os.path.dirname(fimo_file), '{}.sites{}'.format(filename, ext))
    #Split chr-start:end to three columns
    #TODO the _shuf is to handle cases coming from shuffled fasta, this can be generalised and is currently a hack
    split_chr = lambda columnstr: pd.Series(s for s in columnstr.replace('_shuf', '').replace('-', ':').split(':'))
    if fimo_df['sequence name'].str.contains(':|-').sum():
        fimo_df[['chrom', 'chromStart', 'chromEnd']] = fimo_df['sequence name'].apply(split_chr)
        fimo_df.loc[:, ['chromStart', 'chromEnd']] = fimo_df[['chromStart', 'chromEnd']].astype(int)
        fimo_df.loc[:, 'motifStartZeroBased'] = fimo_df.chromStart+fimo_df.start-1
        fimo_df.loc[:, 'motifEndOneBased'] = fimo_df.chromStart+fimo_df.stop
    else:
        fimo_df['chrom'] = fimo_df['sequence name']
        fimo_df['chromStart'] = fimo_df['start']
        fimo_df['chromEnd'] = fimo_df['stop']
        fimo_df.loc[:, 'motifStartZeroBased'] = fimo_df.start-1
        fimo_df['motifEndOneBased'] = fimo_df['stop']


    fimo_df.to_csv(fimo_sites, index=False, sep='\t')
    return fimo_df
