import os
from ..helpers import MocaException
import pandas
from pybedtools import BedTool
import numpy as np
import shutil

__NARROWPEAK_COLUMNS__ = ['chrom', 'chromStart', 'chromEnd',
                         'name', 'score', 'strand',
                         'signalValue', 'p-value', 'q-value']

__BROADPEAK_COLUMNS__ = ['chrom', 'chromStart', 'chromEnd',
                         'name', 'score', 'strand',
                         'signalValue', 'p-value', 'q-value', 'peak']

__BED_TYPES__ = {10: 'narrowPeak',
                 9: 'broadPeak',
                 3: 'questPeak'}

__BED_COLUMN_MAPPING__ = {9: __NARROWPEAK_COLUMNS__,
                          10: __BROADPEAK_COLUMNS__}


class Bedfile(object):
    """Class to crate a bed file object

    Parameters
    ----------
    filepath: string
        Absolute path to bedfile

    genome_table: string
        Absolute path to geonme chromosome size file
    """
    def __init__(self, filepath, genome_table):
        self.filepath = os.path.abspath(filepath)
        self.bed_format = None
        self.extracted_fasta = None
        self.bed = None
        if not os.path.isfile(self.filepath):
            raise MocaException('Bed file {} not found'.format(self.filepath))
        self._read()
        self.bed_format = self.guess_bedformat()
        self.genome_table = os.path.abspath(genome_table)
        assert self.bed_format is not None
        self.determine_peaks()
        self._sort_bed()
        self.write_to_scorefile()

    def _read(self):
        """ Read bedfile as a pandas dataframe
        """
        try:
            self.bed_df = pandas.read_table(self.filepath,
                                        header=None)
        except Exception as e:
            raise MocaException('Error reading bed file {}'.format(self.filepath),
                                'Traceback: {}'.format(e))

    def guess_bedformat(self):
        """Method to guess bed format

        Returns
        -------
        bed_format: string
            BED format

        Example:
            >>> bed_df = Bedfile('file.bed')
            >>> print(bed_df.guess_bed_format())

        """
        self.bed_columns = self.bed_df.columns
        count = len(self.bed_columns)
        try:
            bed_format = __BED_TYPES__[count]
            self.bed_df.columns = __BED_COLUMN_MAPPING__[count]
        except KeyError:
            raise MocaException('Bed file had {} columns. Supported column lengths are {}')
        return bed_format

    def _sort_bed(self):
        """Sort bed by default in descending order of scores
        """
        return self.sort_by(columns=['score'], ascending=False)

    def write_to_scorefile(self):
        """Write bed file as score file

        """
        filename, file_extension = os.path.splitext(self.filepath)
        filename += '.sorted'
        self.scorefile = filename
        self.bed_df.to_csv(self.scorefile, header=False,
                sep='\t',
                index=False,
                columns=['chrom', 'chromStartZeroBased', 'chromEndOneBased', 'name', 'score'])


    def slop_bed(self, flank_length=200):
        """Add flanking sequences to bed file

        Parameters
        ----------
        flank_length: int
            the bed region is expanded in both direction by flank_length number of bases

        Returns
        -------
        slop_bed: dataframe
            Slopped bed data object
        """
        if self.bed is None:
            # Load bed fole into bedtools
            self.bed = BedTool(self.scorefile)
        self.bed.slop(g=self.genome_table, b=flank_length)

    def determine_peaks(self):
        """Converts bed file to a three column file

        Columns: chromosome\t peak_position\t score
        """
        bed_format = self.bed_format
        if bed_format=='narrowPeak':
            """See https://genome.ucsc.edu/FAQ/FAQformat.html#format12
            peak: Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
            chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
            """
            # By default peak position is the  floor of mean of chromStart and chromEnd
            self.bed_df['peak_position'] = np.floor(self.bed_df[['chromStart', 'chromEnd']].mean(axis=1))
            # Peak not called, so chromStart is the peak itself
            self.bed_df.ix[self.bed_df.peak==-1, 'peak_position'] = self.bed_df.ix[self.bed_df.peak==-1, 'chromStart']

            ## Append explicit chromStart and chromEnd position based on peak_position
            self.bed_df['chromStartZeroBased'] = self.bed_df['peak_position'].astype(int)-1
            self.bed_df['chromEndOneBased'] = self.bed_df['peak_position'].astype(int)
        else:
            #TODO
            raise MocaException('Broad region support not implemented')

    def extract_fasta(self, fasta_file, outfile=None):
        """Extract fasta of bed regions

        Parameters
        ----------
        fasta_in: string
            Absolute path to location of reference fasta file

        fasta_out: string
            Path to write extracted fasta sequence

        Returns
        -------
        fasta: string
            Fasta sequence combined
        """
        if self.bed is None:
            # Load bed fole into bedtools
            self.bed = Bedtool(self.scorefile)
        self.extracted_fasta = self.bed.sequence(fi=fasta_file)
        self.temp_fasta = self.extracted_fasta.seqfn
        return self.extracted_fasta

    def sort_by(self, columns=None, ascending=False):
        """Method to sort columns of bedfiles

        Parameters
        ----------
        columns: list
            list of column names to sort by
        ascending: bool
            Sort order(Default: true)

        Returns
        -------
        sorted_bed_df: dataframe
            dataframe with sorted columns
        """
        assert type(columns) is list
        return self.bed_df.sort_values(by=columns, ascending=ascending)
