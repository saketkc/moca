#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""MoCA CLI"""
import os
import re
import sys
from time import sleep
import click
from moca import bedoperations, pipeline
from moca.bedoperations.fimo import get_start_stop_intervals
from moca.helpers import filename_extension
from moca.helpers import read_memefile
from moca.helpers.job_executor import safe_makedir
from moca.plotter import create_plot
from moca import version
from tqdm import tqdm
from click_help_colors import HelpColorsGroup, HelpColorsCommand

conservation_wig_keys = ['phylop', 'gerp', 'phastcons']

class ProgressBar(object):
    def __init__(self, msg_list):
        self.bar = tqdm(msg_list,
                        bar_format='{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]',
                        smoothing=0,
                        miniters=1,
                        mininterval=0.01)

    def show_progress(self, msg):
        self.bar.set_description(msg)
        self.bar.update()
        sleep(0.1)

    def close(self):
        self.bar.close()
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(cls=HelpColorsGroup,
             help_headers_color='yellow',
             help_options_color='green')
@click.version_option(version=version.__version__)
def cli():
    """moca: Motif Conservation Analysis"""
    pass

@cli.command('find_motifs', context_settings=CONTEXT_SETTINGS)
@click.option('--bedfile', '-i',
              help='Bed file input',
              required=True)
@click.option('--oc',
              '-o',
              help='Output Directory',
              required=True)
@click.option('--configuration',
              '-c',
              help='Configuration file',
              required=True)
@click.option('--slop-length',
              default=50,
              help='Flanking sequence length',
              required=True)
@click.option('--flank-motif',
              default=5,
              help='Length of sequence flanking motif',
              required=True)
@click.option('--n-motif',
              default=5,
              help='Number of motifs',
              type=int)
@click.option('--cores',
              '-t',
              default=1,
              help='Number of parallel MEME jobs',
              type=int,
              required=True)
@click.option('--genome-build',
              '-g', '-gb',
              help='Key denoting genome build to use in configuration file',
              required=True)

@click.option('--show-progress',
              help='Print progress',
              is_flag=True)

def find_motifs(bedfile, oc, configuration, slop_length,
                flank_motif, n_motif, cores, genome_build, show_progress):
    """Run meme to locate motifs and create conservation stacked plots"""
    root_dir = os.path.dirname(os.path.abspath(bedfile))
    if not oc:
        moca_out_dir = os.path.join(os.getcwd(), 'moca_output')
    else:
        moca_out_dir = oc
    moca_pipeline = pipeline.Pipeline(configuration)
    genome_data = moca_pipeline.get_genome_data(genome_build)
    genome_fasta = genome_data['fasta']
    genome_table = genome_data['genome_table']

    wigfiles = {}
    for key in list(conservation_wig_keys):
        try:
            wigfiles[key] = genome_data['{}_wig'.format(key)]
        except KeyError:
            pass
    safe_makedir(moca_out_dir)
    bedfile_fn, _ = filename_extension(bedfile)


    if show_progress:
        msg_list = ['Extracting Fasta',
                    'Running MEME',
                    'Running CENTRIMO']
        msg_list_e = ['Generating random Fasta', 'Running fimo random', 'Running fimo main'] + ['Extracting Scores']*len(wigfiles.keys()) + ['Creating PLot']
        msg_list = msg_list + msg_list_e*n_motif
        progress_bar = ProgressBar(msg_list)


    query_train_fasta = os.path.join(moca_out_dir,
                                     bedfile_fn + '_train_flank_{}.fasta'.format(slop_length))
    query_test_fasta = os.path.join(moca_out_dir,
                                    bedfile_fn + '_test_flank_{}.fasta'.format(slop_length))

    if show_progress:
        progress_bar.show_progress('Extracting Fasta')
    bed_o = bedoperations.Bedfile(bedfile, genome_table, moca_out_dir)
    bed_train, bed_test = bed_o.split_train_test_bed(train_peaks_count=500, test_peaks_count=500)

    bed_train_slopped  = bed_o.slop_bed(bed_train, flank_length=slop_length)
    bed_test_slopped  = bed_o.slop_bed(bed_test, flank_length=slop_length)

    bed_o.extract_fasta(bed_train_slopped, fasta_in=genome_fasta, fasta_out=query_train_fasta)
    bed_o.extract_fasta(bed_test_slopped, fasta_in=genome_fasta, fasta_out=query_test_fasta)

    #memechip_out_dir = os.path.join(moca_out_dir, 'memechip_analysis')
    meme_out_dir = os.path.join(moca_out_dir, 'meme_out')
    memechip_out_dir = meme_out_dir
    meme_params = moca_pipeline.get_meme_default_params
    if cores==1:
        re.sub(r' -p*', '', meme_params)
    else:
        re.sub(r'-p*', '-p {}'.format(cores), meme_params)
    if show_progress:
        progress_bar.show_progress('Running MEME')

    #meme_run_out = moca_pipeline.run_memechip(fasta_in=query_fasta, out_dir=memechip_out_dir)
    #
    meme_run_out = moca_pipeline.run_meme(fasta_in=query_train_fasta,
                                          out_dir=meme_out_dir,
                                          strargs=meme_params)
    if meme_run_out['stderr']!='':
        sys.stdout.write('Error running MEME: {}'.format(meme_run_out['stderr']))
        sys.exit(1)
    meme_file = os.path.join(meme_out_dir, 'meme.txt')
    meme_summary = read_memefile(meme_file)

    if show_progress:
        progress_bar.show_progress('Running CENTRIMO')
    centrimo_main_dir = os.path.join(moca_out_dir, 'centrimo_out')
    centrimo_main = moca_pipeline.run_centrimo(meme_file=meme_file,
                                               fasta_in=query_test_fasta,
                                               out_dir=centrimo_main_dir)
    for motif in range(1, meme_summary['total_motifs']+1):
        fimo_rand_dir = os.path.join(memechip_out_dir, 'fimo_random_{}'.format(motif))
        fimo_main_dir = os.path.join(memechip_out_dir, 'fimo_out_{}'.format(motif))

        safe_makedir(fimo_rand_dir)
        random_fasta = os.path.join(fimo_rand_dir, 'random_{}.fa'.format(motif))
        if show_progress:
            progress_bar.show_progress('Generating Random Fasta: {}'.format(motif))
        moca_pipeline.run_fasta_shuffler(fasta_in=query_train_fasta, fasta_out=random_fasta)

        #Random
        if show_progress:
            progress_bar.show_progress('Running FIMO Random')
        fimo_rand = moca_pipeline.run_fimo(motif_file=meme_file,
                                           motif_num=motif,
                                           sequence_file=random_fasta,
                                           out_dir=fimo_rand_dir)
        #Main
        if show_progress:
            progress_bar.show_progress('Running FIMO Main')
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
            if show_progress:
                progress_bar.show_progress('Creating plots')
            sample_score_file = moca_pipeline.save_conservation_scores(main_intervals, wigfile,
                                                                       fimo_main_dir, out_prefix=key)
            control_score_file = moca_pipeline.save_conservation_scores(random_intervals, wigfile,
                                                                        fimo_rand_dir, out_prefix=key)
            sample_score_files.append(sample_score_file)
            control_score_files.append(control_score_file)
        if show_progress:
            progress_bar.show_progress('Creating Plot')
        create_plot(meme_file,
                    bedfile_fn,
                    output_dir=moca_out_dir,
                    centrimo_dir=centrimo_main_dir,
                    motif_number=motif,
                    flank_length=flank_motif,
                    sample_score_files=sample_score_files,
                    control_score_files=control_score_files,
                    reg_plot_titles=[key.capitalize() for key in list(conservation_wig_keys)],
                    annotate=None)

    if show_progress:
        progress_bar.close()


@cli.command('plot', context_settings=CONTEXT_SETTINGS)
@click.option('--meme-dir', '--meme_dir', help='MEME output directory', required=True)
@click.option('--centrimo-dir', '--centrimo_dir', help='Centrimo output directory', required=True)
@click.option('--fimo-dir-sample', '--fimo_dir_sample', help='Sample fimo.txt', required=True)
@click.option('--fimo-dir-control', '--fimo_dir_control', help='Control fimo.txt', required=True)
@click.option('--name', help='Plot title')
@click.option('--flank-motif',
              default=5,
              help='Length of sequence flanking motif',
              required=True)
@click.option('--motif',
              default=1,
              help='Motif number',
              type=int)
@click.option('--oc',
              '-o',
              help='Output Directory',
              required=True)
@click.option('--configuration',
              '-c',
              help='Configuration file',
              required=True)
@click.option('--show-progress',
              help='Print progress',
              is_flag=True)
@click.option('--genome-build',
              '-g', '-gb',
              help='Key denoting genome build to use in configuration file',
              required=True)

def plot(meme_dir, centrimo_dir, fimo_dir_sample, fimo_dir_control, name,
         flank_motif, motif, oc, configuration, show_progress, genome_build):
    """Create stacked conservation plots"""
    if not oc:
        moca_out_dir = os.path.join(os.getcwd(), 'moca_output')
    else:
        moca_out_dir = oc
    if show_progress:
        progress_bar = ProgressBar(['Creating Plots']*len(conservation_wig_keys))
    moca_pipeline = pipeline.Pipeline(configuration)
    meme_file = os.path.join(meme_dir, 'meme.txt')
    genome_data = moca_pipeline.get_genome_data(genome_build)

    wigfiles = {}
    for key in list(conservation_wig_keys):
        try:
            wigfiles[key] = genome_data['{}_wig'.format(key)]
        except KeyError:
            pass

    fimo_control= os.path.join(fimo_dir_control, 'fimo.txt')
    fimo_sample = os.path.join(fimo_dir_sample, 'fimo.txt')
    safe_makedir(moca_out_dir)
    main_intervals = get_start_stop_intervals(fimo_sample, flank_length=flank_motif)
    random_intervals = get_start_stop_intervals(fimo_control, flank_length=flank_motif)

    sample_score_files = []
    control_score_files = []
    for key in list(conservation_wig_keys):
        wigfile = wigfiles[key]
        if show_progress:
            progress_bar.show_progress('Creating plots')
        sample_score_file = moca_pipeline.save_conservation_scores(main_intervals, wigfile,
                                                                    fimo_dir_sample, out_prefix=key)
        control_score_file = moca_pipeline.save_conservation_scores(random_intervals, wigfile,
                                                                    fimo_dir_control, out_prefix=key)
        sample_score_files.append(sample_score_file)
        control_score_files.append(control_score_file)
    create_plot(meme_file,
                name,
                output_dir=moca_out_dir,
                centrimo_dir=centrimo_dir,
                motif_number=motif,
                flank_length=flank_motif,
                sample_score_files=[sample_score_file],
                control_score_files=[control_score_file],
                reg_plot_titles=[key.capitalize() for key in list(conservation_wig_keys)],
                annotate=None)


