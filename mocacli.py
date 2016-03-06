#!/usr/bin/env python
"""MoCA CLI"""

import os
import click
import pyprind
import sys
from moca import bedoperations, wigoperations, pipeline
from moca.bedoperations import fimo_to_sites
from moca.helpers import filename_extension
from moca.helpers import generate_random_fasta
from moca.helpers import read_memefile
from moca.helpers import get_fasta_metadata
from moca.helpers.job_executor import safe_makedir

@click.command()
@click.option('--bedfile', help='Bed file input')
@click.option('--phylop', help='phlyop file input', required=True)
@click.option('--gerp', help='Gerp file input', required=True)
@click.option('--configuration', help='Configuration file', required=True)
@click.option('--genome-table', '-gt', help='Chromosome size file', required=True)
@click.option('--genome-fasta', '-gf', help='Genome fasta', required=True)
@click.option('--flank-length', default=100, help='Flanking sequence length', required=True)

def mocacli(bedfile, genome_table, genome_fasta,
            flank_length, configuration, phylop, gerp):
    """Run moca"""
    root_dir = os.path.dirname(os.path.abspath(bedfile))
    moca_out_dir = os.path.join(root_dir, 'moca_output')
    safe_makedir(moca_out_dir)
    bedfile_fn, _ = filename_extension(bedfile)
    query_fasta = os.path.join(moca_out_dir, bedfile_fn + '_flank_{}.fasta'.format(flank_length))
    bed_df = bedoperations.Bedfile(bedfile, genome_table)
    bed_df.determine_peaks()
    bed_df.slop_bed(flank_length=flank_length)
    bed_df.extract_fasta(fasta_in=genome_fasta, fasta_out=query_fasta)

    moca_pipeline = pipeline.Pipeline(configuration)
    meme_out_dir = os.path.join(moca_out_dir, 'meme_analysis')
    meme_run_out = moca_pipeline.run_meme(fasta_in=query_fasta, out_dir=meme_out_dir)
    meme_file = os.path.join(meme_out_dir, 'meme.txt')

    meme_summary = read_memefile(meme_file)
    fasta_metadata = get_fasta_metadata(query_fasta)
    phylop_wig = wigoperations.WigReader(phylop)
    gerp_wig = wigoperations.WigReader(gerp)
    bar = pyprind.ProgBar(meme_summary['total_motifs']+1, monitor=True, width=180, stream=sys.stdout)

    for motif in range(1, meme_summary['total_motifs']+1):
        fimo_rand_dir = os.path.join(moca_out_dir, 'motif_{}_fimo_analysis_random'.format(motif))
        fimo_main_dir = os.path.join(moca_out_dir, 'motif_{}_fimo_analysis_main'.format(motif))
        safe_makedir(fimo_rand_dir)
        safe_makedir(fimo_main_dir)
        random_fasta = os.path.join(fimo_rand_dir, 'random_{}.fa'.format(motif))
        generate_random_fasta(genome_fasta,
                              genome_table,
                              fasta_metadata['num_seq'],
                              fasta_metadata['len_seq'],
                              random_fasta)

        #Random
        fimo_rand = moca_pipeline.run_fimo(motif_file=meme_file,
                                           motif_num=motif,
                                           sequence_file=random_fasta,
                                           out_dir=fimo_rand_dir)
        #MainMotif
        fimo_main = moca_pipeline.run_fimo(motif_file=meme_file,
                                           motif_num=motif,
                                           sequence_file=query_fasta,
                                           out_dir=fimo_main_dir)

        ##TODO Shoulnd;t this whole thing be a functon????

        fimo_main_file = os.path.join(fimo_main_dir, 'fimo.txt')
        fimo_rand_file = os.path.join(fimo_rand_dir, 'fimo.txt')
        fimo_main_sites = fimo_to_sites(fimo_main_file)
        fimo_rand_sites = fimo_to_sites(fimo_rand_file)

        main_subset = fimo_main_sites[['chrom', 'motifStartZeroBased', 'motifEndOneBased', 'strand']]
        main_intervals = [tuple(x) for x in main_subset.to_records(index=False)]


        rand_subset = fimo_rand_sites[['chrom', 'motifStartZeroBased', 'motifEndOneBased', 'strand']]
        rand_intervals = [tuple(x) for x in rand_subset.to_records(index=False)]


        rand_scores_phylop = phylop_wig.query(rand_intervals)
        rand_scores_gerp = gerp_wig.query(rand_intervals)

        with open(os.path.join(fimo_rand_dir, 'phylop.txt'), 'w') as f:
            f.write('\n'.join(rand_scores_phylop))
        with open(os.path.join(fimo_rand_dir, 'gerp.txt'), 'w') as f:
            f.write('\n'.join(rand_scores_gerp))

        main_scores_phylop = phylop_wig.query(main_intervals)
        main_scores_gerp = gerp_wig.query(main_intervals)

        with open(os.path.join(fimo_main_dir, 'phylop.txt'), 'w') as f:
            f.write('\n'.join(main_scores_phylop))
        with open(os.path.join(fimo_main_dir, 'gerp.txt'), 'w') as f:
            f.write('\n'.join(main_scores_gerp))

        bar.update()


if __name__ == '__main__':
    mocacli()

