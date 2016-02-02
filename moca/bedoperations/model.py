import os
from ..helpers import MocaException
import pandas
from pybedtools import BedTool

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
        self.filepath = filepath
        self.bed_format = None
        if not os.path.isfile(filepath):
            raise MocaException('Bed file {} not found'.format(self.filepath))
        self._read()
        self.bed_format = self.guess_bedformat()
        self.sort_bed()
        self.bed = BedTool(filepath)
        self.genome_table = genome_table
        assert self.bed_Format is not None

    def _read(self):
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
        except KeyError:
            raise MocaException('Bed file had {} columns. Supported column lengths are {}')
        return bed_format

    def slop_bed(self, flank_length=5):
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
        self.bed.slop(g=self.genome_table,
                      b=flank_length
                      )

    def convert_to_scorefile(self):

        """
        filename, file_extension = os.path.splitext(self.filepath)
        filename += '.sorted'
        self.bed_df.to_csv(filename+file_extension,
                           sep='\t',
                           columns=['chrom', 'peak_positions', 'score'],
                           index=False,
                           header=False)
        """
    if filetype=='narrowPeak':
        filter_df1 = df[df.peak.astype(int)==-1]
        filter_df2 = df[df.peak.astype(int)!=-1]
        filter_df1['peak_positions'] = (filter_df1['chromStart'].astype(int)+filter_df1['chromEnd'].astype(int))
        filter_df1['peak_positions'] = [int(x/2) for x in filter_df1['peak_positions'].astype(int)]
        filter_df2['peak_positions'] = filter_df2['chromStart'].astype(int)+filter_df2['peak'].astype(int)
        df = pandas.concat([filter_df1, filter_df2])
    else:
        df['peak_positions'] = (df['chromStart']+df['chromEnd'])
        df['peak_positions'] = [int(x/2) for x in df['peak_positions'].astype(int)]



    def extract_fasta(self, fasta_file):
        """Extract fasta of bed regions
        Parameters
        ----------
        fasta_file: string
            Absolute path to location of fasta file
        Returns
        -------
        fasta: string
            Fasta sequence combined
        """
        self.bed.sequence(fi=fasta_file)


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
        return self.bed_df.sort(columns, ascending)
