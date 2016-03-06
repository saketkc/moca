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
        self.config_file = 'application.cfg'
        self.parsed_config = ConfigurationParser(self.config_file)

    def test_sections(self):
        """Test configuration sections"""
        sections = self.parsed_config.get_all_sections()
        assert 'binaries' in sections

    def test_genomes(self):
        """Test configuration genomes"""
        genomes = self.parsed_config.get_all_genomes()
        assert 'hg19' in genomes


