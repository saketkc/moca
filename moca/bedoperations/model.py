from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object
import os
from moca.helpers import MocaException
from moca.helpers import make_uppercase_fasta
from moca.helpers import path_leaf
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

    output_dir: string
        Output directory
    """
    def __init__(self, filepath, genome_table, output_dir=None):
        self.filepath = os.path.abspath(filepath)
        if not output_dir:
            output_dir = os.path.dirname(self.filepath)
        self.output_dir = output_dir
        self.bed_format = None
        self.extracted_fasta = None
        self.bed = None
        if not os.path.isfile(self.filepath):
            raise MocaException('Bed file {} not found'.format(self.filepath))
        self._read()
        self.bed_format = self.guess_bedformat()
        self.total_peaks = self._count_peaks()
        self.train_peaks_count = None
        self.test_peaks_count = None
        self.genome_table = os.path.abspath(genome_table)
        assert self.bed_format is not None
        self._determine_peaks()
        self._sort_bed()
        self.write_to_scorefile(self.bed_df)

    def _count_peaks(self):
        with open(self.filepath) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

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

    def split_train_test_bed(self, train_peaks_count=500,
                             test_peaks_count=500,
                             train_suffix='train',
                             test_suffix='test'):
        """Split bedfile into ranked training and testing sets

        Training set consists of top `train_peaks_count` peaks
        while testing set will consist of next `test_peaks_count` ranked
        peaks

        Parameters
        ----------
        train_peaks_count: int
            Count of training peaks

        test_peaks_count: int
            Count of testing peaks

        """
        assert train_peaks_count < self.total_peaks
        assert test_peaks_count < self.total_peaks
        self._determine_peaks()
        self._sort_bed()

        self.train_bed_df = self.bed_df[:min(int(train_peaks_count), self.total_peaks)]
        self.test_bed_df = self.bed_df[min(int(test_peaks_count), self.total_peaks):]

        self.train_bed_file = self.write_to_scorefile(self.train_bed_df, out_suffix='train')
        self.test_bed_file = self.write_to_scorefile(self.test_bed_df, out_suffix='test')

        return self.train_bed_file, self.test_bed_file

    def write_to_scorefile(self, bed_df, out_suffix='sorted'):
        """Write bed file as score file

        """
        file_path, file_extension = os.path.splitext(self.filepath)
        filename = path_leaf(file_path)
        filename = os.path.join(self.output_dir, filename+'.{}'.format(out_suffix))
        bed_df.to_csv(filename, header=False,
                      sep='\t',
                      index=False,
                      columns=['chrom', 'peakStartZeroBased',
                               'peakEndOneBased', 'name', 'score'])
        return filename


    def slop_bed(self, bed_in, flank_length=50):
        """Add flanking sequences to bed file

        Parameters
        ----------
        bed_in: string
            Path to input bed file
        flank_length: int
            the bed region is expanded in both direction by flank_length number of bases

        Returns
        -------
        slopped_bed: dataframe
            Slopped bed data object
        """
        inbed = BedTool(bed_in)
        slopped_bed = inbed.slop(g=self.genome_table, b=flank_length)
        filename, extension = os.path.splitext(bed_in)
        filename += '_slop_{}'.format(flank_length) + extension
        slopped_bed.saveas(filename)
        return filename

    def _determine_peaks(self):
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
            raise MocaException('Format should be one of {}'.format(list(__BED_TYPES__.values())))

    def extract_fasta(self, bed_in, fasta_in, fasta_out=None):
        """Extract fasta of bed regions

        Parameters
        ----------
        bed_in: string
            Path to input bed
        fasta_in: string
            Absolute path to location of reference fasta file

        fasta_out: string
            Path to write extracted fasta sequence

        Returns
        -------
        fasta: string
            Fasta sequence combined
        """
        bed = BedTool(bed_in)
        extracted_fasta = bed.sequence(fi=os.path.abspath(fasta_in))
        temp_fasta = extracted_fasta.seqfn
        make_uppercase_fasta(temp_fasta, os.path.abspath(fasta_out))
        os.remove(temp_fasta)

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

    @property
    def get_total_peaks(self):
        return self.total_peaks
