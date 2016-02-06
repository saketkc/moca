"""
test_configurationparser
----------------------------------

Tests for `moca.helpers.configuration_parser` module.
"""

from moca.helpers import ConfigurationParser
import unittest
from moca.helpers import MocaException

class TestConfigParser(unittest.TestCase):

    def setUp(self):
        self.config_file = 'application.cfg'
        self.parsed_config = ConfigurationParser(self.config_file)

    def test_sections(self):
        sections = self.parsed_config.get_all_sections()
        assert 'binaries' in sections

    def test_genomes(self):
        genomes = self.parsed_config.get_all_genomes()
        assert 'hg19' in genomes


