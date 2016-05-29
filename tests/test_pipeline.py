"""
bedoperations_test
----------------------------------

Tests for `moca.bedoperations` module.
"""

import os
import shutil
import unittest
from Bio import SeqIO
from moca.helpers import get_cpu_count
from moca.pipeline import Pipeline
from moca.bedoperations import fimo_to_sites
from moca.helpers import read_memefile

class TestPipeline(unittest.TestCase):
    """Test pipeline """

    def setUp(self):
        """Setup"""
        self.configuration_file = 'tests/data/application.cfg'
        self.meme_fasta = 'tests/data/expected_out/macsPeak.fasta'
        self.pipeline = Pipeline(self.configuration_file)


    def test_meme(self):
        """Test meme runner"""
        if os.path.exists('tests/data/generated_out/meme_analysis'):
            shutil.rmtree('tests/data/generated_out/meme_analysis')
        meme_args = self.pipeline.get_meme_default_params
        output = self.pipeline.run_meme(fasta_in=self.meme_fasta,
                                        out_dir='tests/data/generated_out/meme_analysis',
                                        strargs=meme_args.replace(' -p {}'.format(get_cpu_count()), ''))
        #TODO Check if meme.txt is same and created
        #TODO This check is too stringent, specially if logos are being produced.
        #MEME installation leads to hard coded paths
        print output
        assert output['exitcode'] == 0
        meme_record = read_memefile('tests/data/generated_out/meme_analysis/meme.txt')
        assert meme_record['total_motifs'] == 5

        motif_record1 = meme_record['motif_records'][0]
        motif_record2 = meme_record['motif_records'][1]
        motif_record3 = meme_record['motif_records'][2]

        assert motif_record1.consensus == 'CAGAACGCTGCTGCCAACCCGACCT'
        assert motif_record2.consensus == 'AGCAGA'
        assert motif_record3.consensus == 'CAGTTT'


    def test_fimo(self):
        """Test fimo runner"""
        fimo_output = self.pipeline.run_fimo(motif_file='tests/data/expected_out/meme_analysis/meme.txt',
                                             motif_num=1,
                                             sequence_file=self.meme_fasta,
                                             out_dir='tests/data/generated_out/fimo_analysis',
                                             strargs=None)
        assert fimo_output['exitcode'] == 0

    def test_fimotosites(self):
        """Test fimo_to_sites"""
        fimo_file = 'tests/data/expected_out/fimo_analysis/fimo.txt'
        record_dict = SeqIO.to_dict(SeqIO.parse(open(self.meme_fasta), 'fasta'))
        fimo_df = fimo_to_sites(os.path.abspath(fimo_file))
        for _, row in fimo_df.iterrows():
            record_id = row['sequence name']
            fimo_sequence = row['matched sequence']
            start = row['start']-1 # 0-based start
            end = row['stop'] # 1-based end
            strand = row['strand']
            fasta_sequence = record_dict[record_id].seq
            if strand == '+':
                assert str(fasta_sequence)[start:end] == fimo_sequence
            else:
                assert str(fasta_sequence[start:end].reverse_complement()) == fimo_sequence

    def test_fimostd(self):
        """Test fimo to sites with
        standard sequence names"""
        fimo_file = 'tests/data/fimo.std.txt'
        fimo_df = fimo_to_sites(os.path.abspath(fimo_file))
        for _, row in fimo_df.iterrows():
            record_id = row['sequence name']
            fimo_sequence = row['matched sequence']
            start = row['start']
            end = row['stop'] # 1-based end
            strand = row['strand']
            assert start-1 == row['motifStartZeroBased']
            assert end == row['motifEndOneBased']



    def test_shuffler(self):
        """Smoke test for shuffler"""
        fasta_in = 'tests/data/expected_out/macsPeak.fasta'
        fasta_out = 'tests/data/generated_out/macsPeak.shuffled.fastsa'
        output = self.pipeline.run_fasta_shuffler(fasta_in=fasta_in, fasta_out=fasta_out)
        print output
        with open(fasta_out) as f:
            assert 'chr1' in f.readline()
