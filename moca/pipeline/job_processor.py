"""Job Processor Module
"""
import os
import warnings
from ..helpers import MocaException
from ..helpers import ConfigurationParser
from ..helpers import run_job
from ..helpers import xstr

class Pipeline(object):
    """Generic class to run pipelines

    Parameters
    ----------
    config_file: string
        Optional file input to load all configurations
    """
    def __init__(self, config_file=None):
        self.commands_run = list()
        #TODO This can be removed if config_file is optional
        if not os.path.isfile(xstr(config_file)):
            #TODO This should raise a warning and no xception
            #raise MocaException('Config file {} not found'.format(config_file))
            warnings.warn('No configuration file supplied. Defaults will be used.', UserWarning)
        self.conf = ConfigurationParser(config_file)
        self.meme_default_params = '-dna -mod zoops -nmotifs 3 -minw 6 -maxw 30 -revcomp -nostatus -maxsize 1000000'
        self.meme_strargs = None
        self.meme_location = 'meme'
        self.fimo_default_params = ''
        self.fimo_strargs = None
        self.fimo_location = 'fimo'
        self.shuffler_location = 'fasta-shuffle-letters'
        self.centrimo_args = None
        self.centrimo_location = 'centrimo'
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
            Defualt parametes used: '-dna -revcomp -maxsize 1000000 -nmotifs 3 -p 4'
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
        meme_binary = self.conf.get_binary_path('meme').strip()
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
        fimo_binary = self.conf.get_binary_path('meme')
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
        shuffler_binary = self.conf.get_binary_path('meme')
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

    def run_memechip(self, fasta_in):
        pass

