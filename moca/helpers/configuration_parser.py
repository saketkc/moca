import csv
import ConfigParser
import os

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
        self.allowed_categories = {'genomes': ['wig', 'genome', 'table']}
        self.config_blob = self._read()

    def _read(self):
        #self._check()
        config = ConfigParser.ConfigParser()
        config.read(self.config_file)
        return config

    def _check(self):
        raise NotImplemented

    def get_all_sections(self):
        """Get all sections in config file

        Returns
        -------
        sections: list
            List of all sections
        """
        return self.config_blob.sections()

    def get_all_genomes(self):
        """Get all genome names as specified in config file

        Returns
        -------
        genome_names: list
            List of all genome names
        """
        sections = self.get_all_sections()
        genome_names = filter(lambda x: x.startswith('genome'), sections)
        genome_names = list(map(lambda x: x.replace('genome:',''), genome_names))
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


if __name__ == '__main__':
    config_parser = ConfigurationParser('data/application.cfg')
    print config_parser.get_all_genomes()
