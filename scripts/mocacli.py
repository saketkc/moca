#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""MoCA CLI"""
import os
import numpy as np
import click
from moca import bedoperations, wigoperations, pipeline
from moca.bedoperations import fimo_to_sites
from moca.helpers import filename_extension
from moca.helpers import generate_random_fasta
from moca.helpers import read_memefile
from moca.helpers import get_fasta_metadata
from moca.helpers.job_executor import safe_makedir
from moca.plotter import create_plot
from tqdm import tqdm
from time import sleep
bar = None

def save_score(directory, phylop_wig, gerp_wig, flanking_sites):
    fimo_file = os.path.join(directory, 'fimo.txt')
    print fimo_file
    fimo_sites = fimo_to_sites(fimo_file)

    subset = fimo_sites.loc[:, ['chrom', 'motifStartZeroBased', 'motifEndOneBased', 'strand']]
    subset.loc[:, 'motifStartZeroBased'] = subset['motifStartZeroBased'] - flanking_sites
    subset.loc[:, 'motifEndOneBased'] = subset['motifEndOneBased'] + flanking_sites
    intervals = [tuple(x) for x in subset.to_records(index=False)]

    scores_phylop = phylop_wig.query(intervals)
    scores_gerp = gerp_wig.query(intervals)

    scores_phylop_mean = np.nanmean(scores_phylop, axis=0)
    scores_gerp_mean = np.nanmean(scores_gerp, axis=0)

    np.savetxt(os.path.join(directory, 'phylop.raw.txt'), scores_phylop, fmt='%.4f')
    np.savetxt(os.path.join(directory, 'gerp.raw.txt'), scores_gerp, fmt='%.4f')

    np.savetxt(os.path.join(directory, 'phylop.mean.txt'), scores_phylop_mean, fmt='%.4f')
    np.savetxt(os.path.join(directory, 'gerp.mean.txt'), scores_gerp_mean, fmt='%.4f')

def show_progress(msg):
    bar.set_description(msg)
    bar.update()
    sleep(0.5)

@click.command()
@click.option('--bedfile', help='Bed file input')
@click.option('--phylop', help='phlyop file input', required=True)
@click.option('--gerp', help='Gerp file input', required=True)
@click.option('--configuration', help='Configuration file', required=True)
@click.option('--genome-table', '-gt', help='Chromosome size file', required=True)
@click.option('--genome-fasta', '-gf', help='Genome fasta', required=True)
@click.option('--flank-seq', default=50, help='Flanking sequence length', required=True)
@click.option('--flank-motif', default=5, help='Length of sequence flanking motif', required=True)

def cli(bedfile, genome_table, genome_fasta,
        flank_seq, flank_motif, configuration, phylop, gerp):
    """Run moca"""
    global bar
    root_dir = os.path.dirname(os.path.abspath(bedfile))
    msg_list = ['Extracting Fasta', 'Running meme', 'Reading phylop bigwig', 'Reading gerp bigwig']
    pipeline_msg = ['Generating random fasta', 'Running fimo random', 'Running fimo main', 'Process main scores', 'Process random scores']
    for i in range(0,3):
        msg_list.extend(pipeline_msg)
        bar = tqdm(msg_list, bar_format='{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]',
                   smoothing=0,
                   miniters=1,
                   mininterval=0.01)
    moca_out_dir = os.path.join(root_dir, 'moca_output')
    safe_makedir(moca_out_dir)
    bedfile_fn, _ = filename_extension(bedfile)
    query_fasta = os.path.join(moca_out_dir, bedfile_fn + '_flank_{}.fasta'.format(flank_seq))
    bed_df = bedoperations.Bedfile(bedfile, genome_table)
    bed_df.determine_peaks()
    bed_df.slop_bed(flank_length=flank_seq)
    show_progress('Extracting Fasta')
    bed_df.extract_fasta(fasta_in=genome_fasta, fasta_out=query_fasta)

    moca_pipeline = pipeline.Pipeline(configuration)
    memechip_out_dir = os.path.join(moca_out_dir, 'memechip_analysis')
    meme_out_dir = os.path.join(memechip_out_dir, 'meme_out')
    show_progress('Running meme')

    meme_run_out = moca_pipeline.run_memechip(fasta_in=query_fasta, out_dir=memechip_out_dir)

    meme_file = os.path.join(meme_out_dir, 'meme.txt')
    meme_summary = read_memefile(meme_file)
    fasta_metadata = get_fasta_metadata(query_fasta)

    show_progress('Reading Phylop bigwig')
    phylop_wig = wigoperations.WigReader(phylop)
    show_progress('Reading Gerp bigwig')
    gerp_wig = wigoperations.WigReader(gerp)

    for motif in range(1, meme_summary['total_motifs']+1):
        fimo_rand_dir = os.path.join(memechip_out_dir, 'motif_{}_fimo_analysis_random'.format(motif))
        fimo_main_dir = os.path.join(memechip_out_dir, 'fimo_out_{}'.format(motif))
        safe_makedir(fimo_rand_dir)
        random_fasta = os.path.join(fimo_rand_dir, 'random_{}.fa'.format(motif))
        show_progress('Generating Random Fasta: {}'.format(motif))
        moca_pipeline.run_fasta_shuffler(fasta_in=query_fasta, fasta_out=random_fasta)

        #Random
        show_progress('Running Fimo Random')
        fimo_rand = moca_pipeline.run_fimo(motif_file=meme_file,
                                           motif_num=motif,
                                           sequence_file=random_fasta,
                                           out_dir=fimo_rand_dir)

        show_progress('Processing Scores Random')
        save_score(fimo_rand_dir, phylop_wig, gerp_wig, flank_motif)
        show_progress('Processing Scores Main')
        save_score(fimo_main_dir, phylop_wig, gerp_wig, flank_motif)
        create_plot(meme_file=meme_file,
                    peak_file=os.path.join(root_dir, bedfile_fn + '.sorted'),
                    fimo_file=os.path.join(fimo_main_dir, 'fimo.sites.txt'),
                    sample_phylop_file=os.path.join(fimo_main_dir, 'phylop.mean.txt'),
                    control_phylop_file=os.path.join(fimo_rand_dir, 'phylop.mean.txt'),
                    sample_gerp_file=os.path.join(fimo_main_dir, 'gerp.mean.txt'),
                    control_gerp_file=os.path.join(fimo_rand_dir, 'gerp.mean.txt'),
                    annotate=False,
                    motif_number=motif,
                    flank_length=flank_motif
                    )
    bar.close()

