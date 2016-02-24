import os
from ..helpers import MocaException
from ..helpers import ConfigurationParser
from ..helpers import run_job

class Pipeline(object):
    """Generic class to run pipelines
    Parameters
    ----------
    config_file: string
        Absolute path to configuration file


    """
    def __init__(self, config_file):
        self.commands_run = list()
        if not os.path.isfile(config_file):
            raise MocaException('Config file {} not found'.format(config_file))
        self.conf = ConfigurationParser(config_file)
        self.meme_default_params = '-dna -revcomp -maxsize 1000000 -nmotifs 3'

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
        self.meme_location = self.conf.get_binary_path('meme')
        if not out_dir:
            out_dir = os.path.join(os.path.dirname(fasta_in), 'meme_out')
        stdout, stderr, exitcode = run_job("{} -oc {}".format(self.meme_location, out_dir), cwd=os.path.dirname(out_dir))



    def run_fimo(object):
        pass

    def conservation_analysis(object):
        pass

    def plot_results(object):
        pass
