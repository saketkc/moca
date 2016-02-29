#!/usr/bin/env python
import click
from moca import bedoperations, wigoperations, pipeline
from moca.bedoperations import fimo_to_sites
from moca.helpers import filename_extension
from moca.helpers import generate_random_fasta
from moca.helpers import read_memefile
import os
from moca.helpers.job_executor import safe_makedir

@click.command()
@click.option('--bedfile', help='Bed file input')
@click.option('--configuration', help='Configuration file')
@click.option('--genome-table', '-gt', help='Chromosome size file')
@click.option('--genome-fasta', '-gf', help='Genome fasta')
@click.option('--flank-length', default=100, help='Flanking sequence length')


def mocacli(bedfile, genome_table, genome_fasta, flank_length, configuration):
    root_dir = os.path.dirname(os.path.abspath(bedfile))
    moca_out_dir = os.path.join(root_dir, 'moca_output')
    safe_makedir(moca_out_dir)
    bedfile_fn, ext = filename_extension(bedfile)
    query_fasta = os.path.join(moca_out_dir, bedfile_fn + '_flank_{}.fasta'.format(flank_length))
    bed_df = bedoperations.Bedfile(bedfile, genome_table)
    bed_df.determine_peaks()
    bed_df.slop_bed(flank_length=flank_length)
    bed_df.extract_fasta(fasta_in=genome_fasta, fasta_out=query_fasta)

    moca_pipeline = pipeline.Pipeline(configuration)

    meme_out_dir = os.path.join(moca_out_dir, 'meme_analysis' )
    #def run_meme(self, fasta_in, out_dir=None, strargs=None):
    meme_run_out = moca_pipeline.run_meme(fasta_in=query_fasta, out_dir=meme_out_dir)
    meme_file = os.path.join(meme_out_dir, 'meme.txt')

    meme_summary = read_memefile(meme_file)
    for motif in range(1, meme_summary['num_occurences']+1):
        fimo_rand_dir = os.path.join(moca_out_dir, 'motif_{}_fimo_analysis_random'.format(motif))
        fimo_main_dir = os.path.join(moca_out_dir, 'motif_{}_fimo_analysis_main'.format(motif))
        random_fasta = os.path.join(fimo_random_dir, 'random_{}.fa'.format(motif))
        generate_random_fasta(random_fasta)

        #Random
        fimo_rand = moca_pipeline.run_fimo(motif_file=meme_file, motif_num = motif, sequence_file=random_fasta, out_dir=fimo_rand_dir)
        #MainMotif
        fimo_main = moca_pipeline.run_fimo(motif_file=meme_file, motif_num = motif, sequence_file=query_fasta, out_dir=fimo_main_dir)

        fimo_main_sites = fimo_to_sites(os.path.join(fimo_main_dir,'fimo.txt'))
        fimo_rand_sites = fimo_to_sites(os.path.join(fimo_rand_dir,'fimo.txt'))

        main_subset = fimo_main_sites[['motifStartZeroBased', 'motifEndOneBased', 'strand']]
        main_intervals = [tuple(x) for x in main_subset.values]

        rand_subset = fimo_rand_sites[['motifStartZeroBased', 'motifEndOneBased', 'strand']]
        rand_intervals = [tuple(x) for x in rand_subset.values]






if __name__ == '__main__':
    mocacli()

