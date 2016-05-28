"""
test_configurationparser
----------------------------------

Tests for `moca.helpers.configuration_parser` module.
"""

import unittest
from moca.helpers import ConfigurationParser

class TestConfigParser(unittest.TestCase):
    """Test configurationparser"""

    def setUp(self):
        self.config_file = 'tests/data/application.cfg'
        self.parsed_config = ConfigurationParser(self.config_file)

    def test_sections(self):
        """Test configuration sections"""
        sections = self.parsed_config.get_all_sections()
        assert 'binaries' in sections

    def test_get_section(self):
        """Test configuration sections"""
        section_dict = self.parsed_config.get_section('mongo')
        assert 'mongo_host' in section_dict.keys()
    def test_genomes(self):
        """Test configuration genomes"""
        genomes = self.parsed_config.get_all_genomes()
        assert 'hg19' in genomes


