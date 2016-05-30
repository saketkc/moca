#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""MoCA CLI"""
import os
import click
from moca import bedoperations, pipeline
from moca.bedoperations.fimo import get_start_stop_intervals
from moca.helpers import filename_extension
from moca.helpers import read_memefile
from moca.helpers.job_executor import safe_makedir
from moca.plotter import create_plot
from tqdm import tqdm
from time import sleep
bar = None
conservation_wig_keys = ['phylop', 'gerp', 'phastcons']

def show_progress(msg):
    bar.set_description(msg)
    bar.update()
    sleep(0.5)

@click.command()
@click.option('--bedfile', '-i', help='Bed file input', required=True)
@click.option('--oc', '-o', help='Output Directory')
@click.option('--configuration', '-c', help='Configuration file', required=True)
@click.option('--flank-seq', default=50, help='Flanking sequence length', required=True)
@click.option('--flank-motif', default=5, help='Length of sequence flanking motif', required=True)
@click.option('--n-motif', default=5, help='Number of motifs', type=int)
@click.option('--genome-build', '-g', '-gb',  help='Key denoting genome build to use in configuration file', required=True)


def cli(bedfile, oc, configuration, flank_seq,
        flank_motif, n_motif, genome_build):
    """Run moca"""
    global bar
    root_dir = os.path.dirname(os.path.abspath(bedfile))
    if not oc:
        moca_out_dir = os.path.join(root_dir, 'moca_output')
    else:
        moca_out_dir = oc

    msg_list = ['Extracting Fasta',
                'Running MEME',
                'Running CENTRIMO']
    pipeline_msg = ['Generating random fasta',
                    'Running FIMO on control',
                    'Running FIMO on sample']

    bar = tqdm(msg_list,
               bar_format='{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]',
               smoothing=0,
               miniters=1,
               mininterval=0.01)
    safe_makedir(moca_out_dir)
    bedfile_fn, _ = filename_extension(bedfile)

    moca_pipeline = pipeline.Pipeline(configuration)
    genome_data = moca_pipeline.get_genome_data(genome_build)
    wigfiles = {}
    for key in list(conservation_wig_keys):
        try:
            wigfiles[key] = genome_data['{}_wig'.format(key)]
        except KeyError:
            pass
    genome_fasta = genome_data['fasta']
    genome_table = genome_data['genome_table']

    query_train_fasta = os.path.join(moca_out_dir, bedfile_fn + '_train_flank_{}.fasta'.format(flank_seq))
    query_test_fasta = os.path.join(moca_out_dir, bedfile_fn + '_test_flank_{}.fasta'.format(flank_seq))

    show_progress('Extracting Fasta')
    bed_o = bedoperations.Bedfile(bedfile, genome_table, moca_out_dir)
    bed_train, bed_test = bed_o.split_train_test_bed(train_peaks_count=500, test_peaks_count=500)

    bed_train_slopped  = bed_o.slop_bed(bed_train, flank_length=flank_seq)
    bed_test_slopped  = bed_o.slop_bed(bed_test, flank_length=flank_seq)

    bed_o.extract_fasta(bed_train_slopped, fasta_in=genome_fasta, fasta_out=query_train_fasta)
    bed_o.extract_fasta(bed_test_slopped, fasta_in=genome_fasta, fasta_out=query_test_fasta)

    #memechip_out_dir = os.path.join(moca_out_dir, 'memechip_analysis')
    meme_out_dir = os.path.join(moca_out_dir, 'meme_out')
    memechip_out_dir = meme_out_dir
    show_progress('Running MEME')

    #meme_run_out = moca_pipeline.run_memechip(fasta_in=query_fasta, out_dir=memechip_out_dir)
    #
    meme_run_out = moca_pipeline.run_meme(fasta_in=query_train_fasta, out_dir=meme_out_dir)
    meme_file = os.path.join(meme_out_dir, 'meme.txt')
    meme_summary = read_memefile(meme_file)

    msg_list = []
    for i in range(0, meme_summary['total_motifs']):
        msg_list.extend(pipeline_msg)
        for key in list(conservation_wig_keys):
            msg_list.append('Generating conservation scores')
        msg_list.append('Creating Plots')
    bar = tqdm(msg_list,
               bar_format='{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]',
               smoothing=0,
               miniters=1,
               mininterval=0.01)
    show_progress('Running CENTRIMO')
    centrimo_main_dir = os.path.join(moca_out_dir, 'centrimo_out')
    centrimo_main = moca_pipeline.run_centrimo(meme_file=meme_file,
                                               fasta_in=query_test_fasta,
                                               out_dir=centrimo_main_dir)

    for motif in range(1, meme_summary['total_motifs']+1):
        fimo_rand_dir = os.path.join(memechip_out_dir, 'fimo_random_{}'.format(motif))
        fimo_main_dir = os.path.join(memechip_out_dir, 'fimo_out_{}'.format(motif))

        safe_makedir(fimo_rand_dir)
        random_fasta = os.path.join(fimo_rand_dir, 'random_{}.fa'.format(motif))
        show_progress('Generating Random Fasta: {}'.format(motif))
        moca_pipeline.run_fasta_shuffler(fasta_in=query_train_fasta, fasta_out=random_fasta)

        #Random
        show_progress('Running FIMO Random')
        fimo_rand = moca_pipeline.run_fimo(motif_file=meme_file,
                                           motif_num=motif,
                                           sequence_file=random_fasta,
                                           out_dir=fimo_rand_dir)
        #Main
        show_progress('Running FIMO Main')
        fimo_main = moca_pipeline.run_fimo(motif_file=meme_file,
                                           motif_num=motif,
                                           sequence_file=query_test_fasta,
                                           out_dir=fimo_main_dir)


        fimo_rand_file = os.path.join(fimo_rand_dir, 'fimo.txt')
        fimo_main_file = os.path.join(fimo_main_dir, 'fimo.txt')

        main_intervals = get_start_stop_intervals(fimo_main_file, flank_length=flank_motif)
        random_intervals = get_start_stop_intervals(fimo_rand_file, flank_length=flank_motif)

        sample_score_files = []
        control_score_files = []
        for key in list(conservation_wig_keys):
            wigfile = wigfiles[key]
            show_progress('Creating plots')
            sample_score_file = moca_pipeline.save_conservation_scores(main_intervals, wigfile,
                                                                       fimo_main_dir, out_prefix=key)
            control_score_file = moca_pipeline.save_conservation_scores(random_intervals, wigfile,
                                                                        fimo_rand_dir, out_prefix=key)
            sample_score_files.append(sample_score_file)
            control_score_files.append(control_score_file)

        show_progress('Creating Plot')
        create_plot(meme_file,
                    bedfile_fn,
                    output_dir=os.path.join(oc, 'moca_plots'),
                    centrimo_dir=centrimo_main_dir,
                    motif_number=motif,
                    flank_length=flank_motif,
                    sample_score_files=sample_score_files,
                    control_score_files=control_score_files,
                    reg_plot_titles=[key.capitalize() for key in list(conservation_wig_keys)],
                    annotate=None)

    bar.close()
