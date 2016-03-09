import os
from ..helpers import MocaException
from ..helpers import make_uppercase_fasta
import pandas
from pybedtools import BedTool
import numpy as np

__NARROWPEAK_COLUMNS__ = ['chrom', 'chromStart', 'chromEnd',
                         'name', 'score', 'strand',
                         'signalValue', 'p-value', 'q-value']

__BROADPEAK_COLUMNS__ = ['chrom', 'chromStart', 'chromEnd',
                         'name', 'score', 'strand',
                         'signalValue', 'p-value', 'q-value', 'peak']

__MACSPEAK_COLUMNS__ = __BROADPEAK_COLUMNS__[:5]

__BED_TYPES__ = {10: 'narrowPeak',
                 9: 'broadPeak',
                 5: 'macsPeak'}

__BED_COLUMN_MAPPING__ = {9: __NARROWPEAK_COLUMNS__,
                          10: __BROADPEAK_COLUMNS__,
                          5: __MACSPEAK_COLUMNS__}


class Bedfile(object):
    """Class to crate a bed file object

    Parameters
    ----------
    filepath: string
        Absolute path to bedfile

    genome_table: string
        Absolute path to geonme chromosome size file
    """
    def __init__(self, filepath, genome_table, peaks_to_keep=500):
        self.filepath = os.path.abspath(filepath)
        self.bed_format = None
        self.extracted_fasta = None
        self.bed = None
        self.peaks_to_keep = peaks_to_keep
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
            self.bed_df = pandas.read_table(self.filepath, header=None)
        except Exception as e:
            raise MocaException('Error reading bed file {}\n Traceback: {}'.format(self.filepath, e))

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
            raise MocaException('Bed file had {} columns. Supported column lengths are 5,9 or 10'.format(count))
        return bed_format

    def _sort_bed(self):
        """Sort bed by default in descending order of scores
        """
        return self.sort_by(columns=['score', 'chrom', 'peakStartZeroBased'], ascending=[False, True, True])

    def write_to_scorefile(self):
        """Write bed file as score file

        """
        filename, file_extension = os.path.splitext(self.filepath)
        filename += '.sorted'
        self.scorefile = filename
        self.bed_df[:self.peaks_to_keep].to_csv(self.scorefile, header=False,
                sep='\t',
                index=False,
                columns=['chrom', 'peakStartZeroBased', 'peakEndOneBased', 'name', 'score'])


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
        self.bed = self.bed.slop(g=self.genome_table, b=flank_length)

    def determine_peaks(self):
        """Add extra columns representing peaks

        Parameters
        ----------

        Returns
        --------

        """
        bed_format = self.bed_format
        if bed_format == 'narrowPeak':
            """See https://genome.ucsc.edu/FAQ/FAQformat.html#format12
            peak: Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
            chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
            """
            # By default peka is at 0-based offset from chromStart
            self.bed_df['peak_position'] = self.bed_df['chromStart'] + self.bed_df['peak']
            # Peak not called, so chromStart is the peak itself
            self.bed_df.ix[self.bed_df.peak==-1, 'peak_position'] = self.bed_df.ix[self.bed_df.peak==-1, 'chromStart']

            ## Append explicit chromStart and chromEnd position based on peak_position
            self.bed_df['peakStartZeroBased'] = self.bed_df['peak_position'].astype(int)
            self.bed_df['peakEndOneBased'] = self.bed_df['peak_position'].astype(int)+1
        elif bed_format == 'broadPeak':
            #TODO
            raise MocaException('Broad region support not implemented')
        elif bed_format == 'macsPeak':
            self.bed_df['peak_position'] = np.floor(self.bed_df[['chromStart', 'chromEnd']].mean(axis=1))
            self.bed_df['peakStartZeroBased'] = self.bed_df['peak_position'].astype(int)
            self.bed_df['peakEndOneBased'] = self.bed_df['peak_position'].astype(int)+1
        else:
            raise MocaException('Format should be one of {}'.format(__BED_TYPES__.values()))

    def extract_fasta(self, fasta_in, fasta_out=None):
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
        #if self.bed is None:
            # Load bed fole into bedtools
            #self.bed = BedTool(self.scorefile)
        self.downsampled_bed = self.bed.at(range(0,self.peaks_to_keep))
        self.extracted_fasta = self.downsampled_bed.sequence(fi=os.path.abspath(fasta_in))
        self.temp_fasta = self.extracted_fasta.seqfn
        make_uppercase_fasta(self.temp_fasta, os.path.abspath(fasta_out))
        os.remove(self.temp_fasta)
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
        return self.bed_df.sort_values(by=columns, ascending=ascending, inplace=True)

    def __str__(self):
        return repr(self.bed_df)
