"""Job Processor Module
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import os
import warnings
from ..helpers import ConfigurationParser
from ..helpers import run_job
from ..helpers import xstr
from ..helpers import get_cpu_count
from ..helpers.filename import touch

from ..wigoperations import WigReader
import numpy as np

class Pipeline(ConfigurationParser):
    """Generic class to run pipelines

    Parameters
    ----------
    config_file: string
        Optional file input to load all configurations
    """
    def __init__(self, config_file=None):
        super(Pipeline, self).__init__(config_file)
        self.commands_run = list()
        #TODO This can be removed if config_file is optional
        if not os.path.isfile(xstr(config_file)):
            #TODO This should raise a warning and no xception
            #raise MocaException('Config file {} not found'.format(config_file))
            warnings.warn('No configuration file supplied. Defaults will be used.', UserWarning)
        self.cpu_cores = get_cpu_count()
        self.meme_default_params = '-dna -mod zoops -nmotifs 5 -minw 6 -maxw 30 -revcomp -nostatus -maxsize 1000000 -p {}'.format(self.cpu_cores)
        self.meme_strargs = None
        self.meme_location = 'meme'
        self.fimo_default_params = ''
        self.fimo_strargs = None
        self.fimo_location = 'fimo'
        self.shuffler_location = 'fasta-shuffle-letters'
        self.centrimo_args = None
        self.centrimo_location = 'centrimo'
        self.memechip_default_params = '-dna -meme-mod zoops -meme-nmotifs 5 -meme-minw 6 -meme-maxw 30 -meme-maxsize 1000000 -meme-p {}'.format(self.cpu_cores)
        self.memechip_args = None
        self.memechip_location = 'meme-chip'
        self.commands_run = []

    def run_meme(self, fasta_in, out_dir=None, strargs=None):
        """Run meme
        Run meme on a given input fasta

        Parameters
        ---------
        strargs: string
            A concatenated string containing parameters as would be passed to standalone meme
            Defualt parametes used: '-dna -revcomp -maxsize 1000000 -nmotifs 5 -p 4'
            To modify: Pipeline.meme_default_params
        fasta_in: string
            Location of the fasta file
        out_dir: string
            Location to write all meme analysis output

        Returns
        -------
        meme_out: string
            Location of meme output
        """
        self.meme_strargs = strargs
        if not self.meme_strargs:
            self.meme_strargs = self.meme_default_params
        meme_binary = self.get_binary_path('meme').strip()
        if not meme_binary or meme_binary == '':
            # Use meme from envirnonment
            meme_binary = 'meme'
        else:
            #  Use absolute path meme
            meme_binary += '/meme'
        self.meme_location = meme_binary
        if not out_dir:
            out_dir = os.path.join(os.path.dirname(fasta_in), 'meme_out')
        out_dir = os.path.abspath(out_dir)
        cmd = '{} {} -oc {} {}'.format(self.meme_location, self.meme_strargs,
                                       out_dir, os.path.abspath(fasta_in))
        stdout, stderr, exitcode = run_job(cmd=cmd,
                                           cwd=os.path.dirname(out_dir))

        output = {'out_dir': out_dir, 'stdout': stdout,
                  'stderr': stderr, 'exitcode': exitcode,
                  'cmd': cmd}
        self.commands_run.append({'cmd': cmd, 'metadata': output})
        return output

    def run_fimo(self, motif_file, motif_num, sequence_file, out_dir=None, strargs=None):
        """Run fimo to find out locations where motif occurs

        Parameters
        ---------
        motif_file: str
            Path to meme.txt
        motif_num: int
            Motif number to investigate
        sequence_file: str
            Path to sequence file
        out_dir: str
            Location to output fimo results
        strargs: str
            string arguments as would be passed to fimo commandline
        """
        #TODO This code is same as in the above fmethod. Make this a separate method?
        self.fimo_strargs = strargs
        if not out_dir:
            out_dir = os.path.join(os.path.dirname(motif_file), 'fimo_out')
        self.fimo_strargs = xstr(self.fimo_strargs) + ' --motif {} -oc {}'.format(motif_num, os.path.abspath(out_dir))
        fimo_binary = self.get_binary_path('meme')
        if not fimo_binary or fimo_binary == '':
            # Use meme from envirnonment
            fimo_binary = 'fimo'
        else:
            #  Use absolute path meme
            fimo_binary += '/fimo'
        self.fimo_location = fimo_binary
        cmd = '{}{} {} {}'.format(self.fimo_location, self.fimo_strargs,
                                  os.path.abspath(motif_file),
                                  os.path.abspath(sequence_file))
        stdout, stderr, exitcode = run_job(cmd=cmd,
                                           cwd=os.path.dirname(out_dir))
        output = {'out_dir': out_dir, 'stdout': stdout,
                  'stderr': stderr, 'exitcode': exitcode,
                  'cmd': cmd}
        self.commands_run.append({'cmd': cmd, 'metadata': output})
        return output

    def run_fasta_shuffler(self, fasta_in, fasta_out):
        """Run fasta-dinucleotide-shuffle to generate random fasta"""
        shuffler_binary = self.get_binary_path('meme')
        if not shuffler_binary or shuffler_binary == '':
            # Use meme from envirnonment
            shuffler_binary = 'fasta-shuffle-letters'
        else:
            #  Use absolute path meme
            shuffler_binary += '/fasta-shuffle-letters'
        self.shuffler_location = shuffler_binary
        cmd = '{} -kmer 2 -dna {} '.format(self.shuffler_location,
                                             os.path.abspath(fasta_in))
        stdout, stderr, exitcode = run_job(cmd=cmd,
                                           cwd=os.path.dirname(fasta_out))
        with open(os.path.abspath(fasta_out), 'w') as f:
            f.write(stdout)
        output = {'stdout': stdout,
                  'stderr': stderr,
                  'exitcode': exitcode,
                  'cmd': cmd}
        return output

    def run_centrimo(self, fasta_in, meme_file, out_dir):
        """Run centrimo

        Parameters
        ----------
        fasta_in: string
            Path to input fasta
        meme_file: string
            Path to MEME's motif file
        out_dir: string
            Output directory
        Returns
        -------
        output: dict
            A dictionary with 'stderr,stdout,cmd,exitcode,out_dir'

        """
        centrimo_binary = self.get_binary_path('meme').strip()
        if not centrimo_binary or centrimo_binary=='':
            centrimo_binary = 'centrimo'
        else:
            centrimo_binary += '/centrimo'
        self.centrimo_location = centrimo_binary
        if not out_dir:
            out_dir = os.path.join(os.path.abspath(os.path.join(os.path.dirname(fasta_in), os.pardir)), 'centrimo_out')
        else:
            out_dir = os.path.abspath(out_dir)

        cmd = '{} -oc {} {} {}'.format(self.centrimo_location, out_dir, os.path.abspath(fasta_in), os.path.abspath(meme_file))
        stdout, stderr, exitcode = run_job(cmd=cmd,
                                           cwd=os.path.dirname(out_dir))

        output = {'out_dir': out_dir, 'stdout': stdout,
                  'stderr': stderr, 'exitcode': exitcode,
                  'cmd': cmd}
        self.commands_run.append({'cmd': cmd, 'metadata': output})
        return output


    def run_memechip(self, fasta_in, out_dir=None, strargs=None):
        """Run meme-chip
        Run meme-chip on a given input fasta

        Parameters
        ---------
        strargs: string
            A concatenated string containing parameters as would be passed to standalone meme
            Defualt parametes used: '-dna -revcomp -maxsize 1000000 -nmotifs 5 -p 4'
            To modify: Pipeline.meme_default_params
        fasta_in: string
            Location of the fasta file
        out_dir: string
            Location to write all meme analysis output

        Returns
        -------
        output: dict
            A dictionary with 'stderr,stdout,cmd,exitcode,out_dir'
        """
        self.memechip_strargs = strargs
        if not self.memechip_strargs:
            self.memechip_strargs = self.memechip_default_params
        meme_binary = self.get_binary_path('meme').strip()
        if not meme_binary or meme_binary == '':
            # Use meme from envirnonment
            meme_binary = 'meme-chip'
        else:
            #  Use absolute path meme
            meme_binary += '/meme-chip'
        self.meme_location = meme_binary
        if not out_dir:
            out_dir = os.path.join(os.path.dirname(fasta_in), 'memechip_out')
        out_dir = os.path.abspath(out_dir)
        cmd = '{} {} -oc {} {}'.format(self.meme_location, self.memechip_strargs,
                                       out_dir, os.path.abspath(fasta_in))
        stdout, stderr, exitcode = run_job(cmd=cmd,
                                           cwd=os.path.dirname(out_dir))

        output = {'out_dir': out_dir, 'stdout': stdout,
                  'stderr': stderr, 'exitcode': exitcode,
                  'cmd': cmd}
        self.commands_run.append({'cmd': cmd, 'metadata': output})
        return output

    @staticmethod
    def save_conservation_scores(intervals, wig_file, out_directory, out_prefix='phylop'):
        """Extract and save conservation scores
        Parameters
        ----------
        intervals: list
            list of tuple

        wig_file: string
            Wig file location

        out_directory: string
            Out file directory

        out_prefix: string
            Out file prefix
        """
        wig = WigReader(wig_file)
        conservation_scores = wig.query(intervals)
        if np.any(conservation_scores):
            conservation_scores_mean = np.nanmean(conservation_scores, axis=0)
            np.savetxt(os.path.join(out_directory, '{}.raw.txt'.format(out_prefix)),
                       conservation_scores, fmt='%.4f')
            np.savetxt(os.path.join(out_directory, '{}.mean.txt'.format(out_prefix)),
                       conservation_scores_mean, fmt='%.4f')
        else:
            touch(os.path.join(out_directory, '{}.mean.txt'.format(out_prefix)))
            touch(os.path.join(out_directory, '{}.raw.txt'.format(out_prefix)))
        return os.path.join(out_directory, '{}.mean.txt'.format(out_prefix))

    @property
    def get_meme_default_params(self):
        return self.meme_default_params

    @property
    def get_memechip_default_params(self):
        return self.memechip_default_params

    @property
    def get_fimo_default_params(self):
        return self.fimo_default_params


