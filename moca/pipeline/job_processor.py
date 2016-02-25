"""Job Processor Module
"""
import os
from ..helpers import MocaException
from ..helpers import ConfigurationParser
from ..helpers import run_job

xstr = lambda s: s or ""

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
        if not os.path.isfile(config_file):
            raise MocaException('Config file {} not found'.format(config_file))
        self.conf = ConfigurationParser(config_file)
        self.meme_default_params = '-dna -revcomp -maxsize 1000000 -nmotifs 3'
        self.meme_strargs = None
        self.meme_location = 'meme'
        self.commands_run = []

    def run_meme(self, fasta_in, out_dir=None, strargs=None):
        """Run meme
        Arguments
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
        meme_binary = self.conf.get_binary_path('meme')
        if not meme_binary:
            # Use meme from envirnonment
            meme_binary = 'meme'
        else:
            #  Use absolute path meme
            meme_binary += '/meme'
        self.meme_location = meme_binary
        if not out_dir:
            out_dir = os.path.join(os.path.dirname(fasta_in), 'meme_out')
        out_dir = os.path.abspath(out_dir)
        cmd = "{} {} -oc {} {}".format(self.meme_location, self.meme_strargs, out_dir, os.path.abspath(fasta_in))
        stdout, stderr, exitcode = run_job(cmd=cmd,
                                           cwd=os.path.dirname(out_dir))

        output = {'out_dir': out_dir, 'stdout': stdout, 'stderr': stderr, 'exitcode': exitcode, 'cmd': cmd}
        self.commands_run.append({'cmd': cmd, 'metadata': output})
        return output
