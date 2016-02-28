"""Convert fimo.txt to a bed file"""
import pandas as pd
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
    fimo_df = pd.read_table(fimo_file)
    #Split chr-start:end to three columns
    split_chr = lambda columnstr: pd.Series(s for s in columnstr.replace('-', ':').split(':'))
    fimo_df[['chrom', 'chromStart', 'chromEnd']] = fimo_df['sequence name'].apply(split_chr)
    fimo_df[['chromStart', 'chromEnd']] = fimo_df[['chromStart', 'chromEnd']].astype(int)
    fimo_df.motifStartZeroBased = fimo_df.chromStart+fimo_df.start-1
    fimo_df.motifEndOneBased = fimo_df.chromEnd+fimo_df.stop
    # Seq((str(fasta[3].seq)[16:40])).reverse_complement()
    return fimo_df





