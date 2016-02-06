import subprocess
from ..helpers import MocaException
from ..helpers import ConfigurationParser

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

    def run_meme(object):
        pass

    def run_fimo(object):
        pass

    def conservation_analysis(object):
        pass

    def plot_results(object):
        pass
