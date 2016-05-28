from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import object
import csv
from ..helpers import MocaException
import configparser
import os

def is_executable(fpath):
    """Check if binary is executable

    Source: http://stackoverflow.com/a/377028/756986
    """
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


class ConfigurationParser(object):
    """
    Parse config file

    Parameters
    ----------
    config_file: str
        absolute path to configuration file

    """

    def __init__(self, config_file):
        self.config_file = config_file
        if not os.path.isfile(config_file):
            raise MocaException('Config file {} not found'.format(config_file))
        self.allowed_categories = {'genomes': ['wig', 'genome', 'table']}
        self.config_blob = self._read()

    def _read(self):
        #self._check()
        config = configparser.ConfigParser()
        config.read(self.config_file)
        return config

    def _check(self):
        raise NotImplementedError

    def get_all_sections(self):
        """Get all sections in config file

        Returns
        -------
        sections: list
            List of all sections
        """
        return self.config_blob.sections()

    def get_section(self, section_name):
        """Get all keys inside a particular section

        Parameters
        ----------
        section_name: str
            Section Name

        Returns
        -------
        section_dict: dict
            Dict of all keys inside header
        """
        section_dict = dict(self.config_blob.items(section_name))
        return section_dict

    def get_all_genomes(self):
        """Get all genome names as specified in config file

        Returns
        -------
        genome_names: list
            List of all genome names
        """
        sections = self.get_all_sections()
        genome_names = [x for x in sections if x.startswith('genome')]
        genome_names = list([x.replace('genome:','') for x in genome_names])
        return genome_names

    def get_genome_data(self, genome_name):
        """ Return all files associated with a genome in the  config file

        Parameters
        ----------
        genome_name: str
            Genome name as specified in 'genome:<genome_name>' format

        Returns
        -------
        genome_dict: str
            Dictionary specifying the files associated with a genome
            Format: {'genome_name': {'<wig_prefix>': wig_prefix_path, 'genome_table': genome_table_path}}
        """
        genome_dict = dict(self.config_blob.items('genome:{}'.format(genome_name)))
        return genome_dict

    def get_binary_path(self, binary_name):
        """ Returns absolute path to installed binaries

        Parameters
        ----------
        binary_name: string
            Program name

        Returns
        -------
        binary_path: string
            Absolute path to installed binary
        """
        binary_path = dict(self.config_blob.items('binaries'))[binary_name]
        if binary_name!='meme':
            #Since meme invovles a lot of subprograms,
            #we avoid checking executalble status
            if not is_executable(binary_path):
                raise MocaException('{} is not executable.'.format(binary_path))
        else:
            if not is_executable(os.path.join(binary_path,'meme')):
                pass
                #TODO This should also raise warning
                #raise MocaException('{}/meme is not executable'.format(binary_path))
        return binary_path

