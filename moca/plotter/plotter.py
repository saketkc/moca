#!/usr/bin/env python
"""
Generate conservation plots
"""
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='monospace')
plt.rcParams.update({'axes.titlesize': 'small'})
plt.rcParams.update({'backend' : 'Agg'})
import csv
import sys
import numpy as np
from scipy import stats
from scipy.stats.stats import pearsonr
import statsmodels.api as sm
from Bio import motifs
import ntpath
import os
from math import log
import matplotlib.gridspec as gridspec
from pylab import setp
from scipy.stats import binom
from matplotlib.font_manager import FontProperties
flankingstemcolor='y'
import json
import click

ENRICHMENT_SEQ_LENGTH = 401

a=39.33333333
__scale__ = 1#0.51#2.5
# Use 'pssm' or 'counts' to score:
score_type = 'counts' #'pssm'
bases = ['A', 'T', 'G', 'C']
histogram_nbins = 20

## Plot parameters
linewidth=3
plot_linewidth=5
markeredgewidth=3.5
greycolor='0.65'
markersize=13
fontsize=20
legend_fontsize = 26
pointsize=20
dpi=300
tickpad = 20
ticklength = 10


plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

##  Text position
txtx=0.02
txty=0.09

legend_xmultiplier = 1.4
legend_ymultiplier = 1


def path_leaf(path):
    """
    Returns the parent path, and the childname
    given an absolute path.
    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def save_figure(fig, ax, name, form):
    """
    Save a given subplot(ax)
    """
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig('{}.{}'.format(name, form), format=form, dpi=dpi, bbox_inches=extent.expanded(1.2, 1.1))
    if form!='png':
        fig.savefig('{}.{}'.format(name, 'png'), format='png', dpi=dpi, bbox_inches=extent.expanded(1.2, 1.1))

def position_wise_profile(counts_dict, length):
    """
    Convert base to position wise profile
    """
    profile = map(dict, zip(*[[(k, v) for v in value] for k, value in counts_dict.items()]))
    return profile

def find_max_occurence(profile, max_count=2):
    """
    Return profile with base corresponding to max scores[CHECK!]
    """
    sorted_profile = []
    for p in profile:
        sorted_profile.append(sorted(p.items(), key=lambda x:x[1]))
    for i,p in enumerate(sorted_profile):
        sorted_profile[i] = p[-max_count:]
    return sorted_profile

def perform_t_test(a,b):
    return stats.ttest_ind(a,b)

def read_peaks(peak_file):
    """
    Peak file reader
    """
    peaks_by_chr = {}
    scores_by_chr = {}
    with open(peak_file,'r') as f:
        for line in f:
            split = line.split('\t')
            chromosome = split[0]
            position = int(split[1])
            score = float(split[4])
            if chromosome not in peaks_by_chr.keys():
                peaks_by_chr[chromosome] = []
                scores_by_chr[chromosome] = []
            peaks_by_chr[chromosome].append(position)
            scores_by_chr[chromosome].append(score)
    return peaks_by_chr, scores_by_chr

def read_fimo_file(fimo_file):
    """
    Fimo file reader
    """
    f = open(fimo_file, 'r')
    reader = csv.reader(f, delimiter='\t')
    reader.next()
    peaks_by_chr = {}
    strand_by_chr = {}
    for row in reader:
        chromosome = row[-3]
        start = int(row[-2])
        end = int(row[-1])
        motif_center = (start+end)/2.0
        strand = row[4]

        if chromosome not in peaks_by_chr:
            peaks_by_chr[chromosome] = []
            strand_by_chr[chromosome] = []
        peaks_by_chr[chromosome].append(motif_center)
        strand_by_chr[chromosome].append(strand)

    return peaks_by_chr, strand_by_chr

def get_motif_distances(peak_file, fimo_file):
    """
    Get distance between motif center and peak position
    """
    peaks_by_chr, scores_by_chr = read_peaks(peak_file)
    fimo_data, strand_by_chr = read_fimo_file(fimo_file)
    chr_wise_distances = {}
    for chr in peaks_by_chr.keys():
        peak_positions = peaks_by_chr[chr]
        try:
            fimo_positions = fimo_data[chr]
        except KeyError:
            ## There is no motif site correpsonding to the chromosome of the peak
            continue
        distances = []
        ## For each peak position we calculate the distance
        ## from its centere to all known  fimo sites in the same chromsome
        for peak_pos in peak_positions:
            for index, fimo_pos in enumerate(fimo_positions):
                strand = strand_by_chr[chr][index]
                if strand=='+':
                    distances.append(fimo_pos-peak_pos)
                else:
                    distances.append(-(fimo_pos-peak_pos))
        chr_wise_distances[chr] = distances
    ## Flatten list of lists
    all_distances = [item for sublist in chr_wise_distances.values() for item in sublist]
    return all_distances

def create_plot(meme_file, motif_number, flanking_sites, sample_phylop_file, control_phylop_file, sample_gerp_file, control_gerp_file, peak_file, fimo_file, annotate):
    handle = open(meme_file)
    records = motifs.parse(handle, 'meme')
    record = records[motif_number-1]
    num_occurrences = getattr(record, 'num_occurrences', 'Unknown')
    sample_phylo_data = None
    control_phylo_data = None
    sample_gerp_data = None
    control_gerp_data = None
    annotate_dict = None
    if annotate == "" or annotate == ' ':
        annotate = None
    elif annotate:
        with open(annotate) as f:
            annotate_dict = json.load(f)


    handle =  open(sample_phylop_file, 'r')
    sample_phylo_data = csv.reader(handle, delimiter='\t')

    handle = open(control_phylop_file, 'r')
    control_phylo_data = csv.reader(handle, delimiter='\t')


    if sample_gerp_file and control_gerp_file:

        handle = open(sample_gerp_file, 'r')
        sample_gerp_data = csv.reader(handle, delimiter='\t')

        handle = open(control_gerp_file, 'r')
        control_gerp_data = csv.reader(handle, delimiter='\t')

    sample_phylo_scores = []
    for line in sample_phylo_data:
        sample_phylo_scores.append(float(line[0]))
    control_phylo_scores = []
    for line in control_phylo_data:
        control_phylo_scores.append(float(line[0]))

    print('LENFTH sample_phylop: {}'.format(len(sample_phylo_scores)))
    if sample_gerp_data:
        sample_gerp_scores = []
        for line in sample_gerp_data:
            sample_gerp_scores.append(float(line[0]))
        control_gerp_scores = []
        for line in control_gerp_data:
            control_gerp_scores.append(float(line[0]))

    assert (len(sample_phylo_scores) == len(control_phylo_scores))

    handle.close()
    profile = position_wise_profile(getattr(record, score_type), record.length)
    max_occur = find_max_occurence(profile, max_count=1)
    ## motif_scores is tn array of scores of the max  occuring base at each position of the motif
    motif_scores = []
    for position in max_occur:
        motif_scores.append(position[0][1])

    motif_scores = np.asarray(motif_scores)
    sample_phylo_scores = np.asarray(sample_phylo_scores)
    control_phylo_scores = np.asarray(control_phylo_scores)
    if sample_gerp_data:
        sample_gerp_scores = np.asarray(sample_gerp_scores)
        control_gerp_scores = np.asarray(control_gerp_scores)

    motif_junk = [0 for i in range(0, flanking_sites)]
    motif_junk = np.asarray(motif_junk)
    motif_concat = np.concatenate((motif_junk, motif_scores))
    motif_concat = np.concatenate((motif_concat, motif_junk))


    ##Mean of flanking sites
    ms_p = np.mean(np.concatenate((sample_phylo_scores[0:flanking_sites], sample_phylo_scores[-flanking_sites:])))
    mc_p = np.mean(np.concatenate((control_phylo_scores[0:flanking_sites], control_phylo_scores[-flanking_sites:])))


    if sample_gerp_data:
        ms_g = np.mean(np.concatenate((sample_gerp_scores[0:flanking_sites], sample_gerp_scores[-flanking_sites:])))
        mc_g = np.mean(np.concatenate((control_gerp_scores[0:flanking_sites], control_gerp_scores[-flanking_sites:])))
        flanking_sample_gerp_scores = np.concatenate((sample_gerp_scores[0:flanking_sites], sample_gerp_scores[-flanking_sites:]))
        flanking_control_gerp_scores = np.concatenate((control_gerp_scores[0:flanking_sites], control_gerp_scores[-flanking_sites:]))
        motif_control_gerp_scores = control_gerp_scores[flanking_sites:-flanking_sites]
        motif_sample_gerp_scores = sample_gerp_scores[flanking_sites:-flanking_sites]

    flanking_sample_phylo_scores = np.concatenate((sample_phylo_scores[0:flanking_sites], sample_phylo_scores[-flanking_sites:]))
    flanking_control_phylo_scores = np.concatenate((control_phylo_scores[0:flanking_sites], control_phylo_scores[-flanking_sites:]))
    motif_control_phylo_scores = control_phylo_scores[flanking_sites:-flanking_sites]
    motif_sample_phylo_scores = sample_phylo_scores[flanking_sites:-flanking_sites]

    if flanking_sites>0:
        shifted_sample_phylo_scores = sample_phylo_scores[flanking_sites:-flanking_sites]-ms_p
        shifted_control_phylo_scores = control_phylo_scores[flanking_sites:-flanking_sites]-mc_p
        if sample_gerp_data:
            shifted_sample_gerp_scores = sample_gerp_scores[flanking_sites:-flanking_sites]-ms_g
            shifted_control_gerp_scores = control_gerp_scores[flanking_sites:-flanking_sites]-mc_g
    else:
        shifted_sample_phylo_scores = sample_phylo_scores
        shifted_control_phylo_scores = control_phylo_scores
        if sample_gerp_data:
            shifted_sample_gerp_scores = sample_gerp_scores
            shifted_control_gerp_scores = control_gerp_scores

    pr_p = pearsonr(motif_scores, motif_sample_phylo_scores)
    if sample_gerp_data:
        pr_g = pearsonr(motif_scores, motif_sample_gerp_scores)

    ## H_0: Mean phylop scores for motif sites and flanking sites are the same
    ## H_!: Mean phylop score for motif sites > Mean phylop score of flanking sites
    ## NOTE: the perform_t_test functions returns a 2 tailed p-value forn independet t-test with unequal sample size, eqaul variances

    T_deltaphylop, p_deltaphylop = perform_t_test(motif_sample_phylo_scores, flanking_sample_phylo_scores)
    delta_phylop = np.mean(motif_sample_phylo_scores)-np.mean(flanking_sample_phylo_scores)#-shifted_control_phylo_scores)
    if sample_gerp_data:
        T_deltagerp, p_deltagerp = perform_t_test(motif_sample_gerp_scores, flanking_sample_gerp_scores)
        delta_gerp = np.mean(motif_sample_gerp_scores)-np.mean(flanking_sample_gerp_scores)
        if T_deltagerp<0:
            p_deltagerp = 1-p_deltagerp*0.5
        else:
            p_deltagerp = p_deltagerp*0.5


    if T_deltaphylop<0:
        p_deltaphylop = 1-p_deltaphylop*0.5
    else:
        p_deltaphylop = p_deltaphylop*0.5


    ## Ordinary least square fit for phylop scores and motif_scores
    reg_phylop_sample = sm.OLS(motif_sample_phylo_scores,sm.add_constant(motif_scores)).fit()
    if (len(reg_phylop_sample.params)<2):
        y_reg_phylop_sample = motif_scores
    else:
        y_reg_phylop_sample = motif_scores*reg_phylop_sample.params[1]+reg_phylop_sample.params[0]
    reg_phylop_control = sm.OLS(motif_control_phylo_scores,sm.add_constant(motif_scores)).fit()
    if (len(reg_phylop_control.params)<2):
        y_reg_phylop_control = motif_scores
    else:
        y_reg_phylop_control = motif_scores*reg_phylop_control.params[1]+reg_phylop_control.params[0]


    if sample_gerp_data:
        reg_gerp_sample = sm.OLS(motif_sample_gerp_scores,sm.add_constant(motif_scores)).fit()
        if (len(reg_gerp_sample.params)==1):
            y_reg_gerp_sample = motif_scores
        else:
            y_reg_gerp_sample = motif_scores*reg_gerp_sample.params[1]+reg_gerp_sample.params[0]

        reg_gerp_control = sm.OLS(motif_control_gerp_scores,sm.add_constant(motif_scores)).fit()
        if (len(reg_gerp_control.params)==1):
            y_reg_gerp_control = motif_scores
        else:
            y_reg_gerp_control = motif_scores*reg_gerp_control.params[1]+reg_gerp_control.params[0]

    motif = record
    motif_length = motif.length
    meme_dir = os.path.dirname(meme_file)
    X = [40+15] ## this is by trial and error, the position for the first base logo
    logo = plt.imread(os.path.join(meme_dir, 'logo{}.png'.format(motif_number)))
    ## Generate all other X coordinates
    fs = flanking_sites
    for j in range(1,len(motif)+2*fs):
        t = X[j-1]+a+1.9
        X.append(t)
    motif_bits = []
    for i in range(0, motif.length):
        s = 0
        for base in bases:
            s = s + -motif.pwm[base][i]*log(motif.pwm[base][i], 2) if motif.pwm[base][i]!=0 else s
            s = 2-s
        motif_bits.append(s)

    y_phylop_pixels = [__scale__*x for x in sample_phylo_scores]#[fs:-fs]]#[flanking_sites:-flanking_sites]]


    ##FIXME This is a big dirty hacl to get thegenerate plots for the Reverse complement logo too
    logo_name =['logo{}.png'.format(motif_number), 'logo_rc{}.png'.format(motif_number)]
    for ln in logo_name:
        if 'rc'in ln:
            y_phylop_pixels.reverse()
        logo = plt.imread(os.path.join(meme_dir, ln))
        height_px = logo.shape[0] # Should be 212

        if sample_gerp_data:
            if annotate:
                total_px = X[-1]+8*height_px+140
                right = (8*height_px+10+140-0.2*height_px)/total_px
            else:
                total_px = X[-1]+6*height_px+140
                right = (6*height_px+10+140-0.2*height_px)/total_px
        else:
            if annotate:
                total_px = X[-1]+6*height_px+140
                right = (6*height_px+10+140-0.2*height_px)/total_px
            else:
                total_px = X[-1]+4*height_px+140
                right = (4*height_px+10+140-0.2*height_px)/total_px

        figsize=(total_px/100,(2*height_px)/100+0.6)

        gs =  gridspec.GridSpec(2, 1)#, width_ratios=[1, right], height_ratios=[1,1])
        gs.update(top=1.0, bottom=0.14, left=0.08, right=1-right)#, right=0.8)#, left=0.06)#, right=right, wspace=0.025, hspace=0.03, wd)
        f = plt.figure(figsize=figsize, dpi=dpi, facecolor='w', edgecolor='k')

        # ax => Logo
        # stem_plot => Trend
        # gerp_scatter_plot => Phylop
        # enrichment_plot => Gerp
        logo_plot = plt.Subplot(f, gs[0])
        ##TODO Check this
        if motif_length>45:
            XSCALE_FACTOR = motif_length/1.9
            z=2
        elif motif_length>40:
            XSCALE_FACTOR = motif_length/2.25
            z=2.5
        elif motif_length>36:
            XSCALE_FACTOR = motif_length/1.95
            z=2
        elif motif_length>21:
            XSCALE_FACTOR = motif_length/5
            z=3
        else:
            XSCALE_FACTOR = 4.5
            z=3

        logo_plot.imshow(logo, extent=[40+15+z*(a+1.9),logo.shape[1]+15+XSCALE_FACTOR*(a+1.9),0,logo.shape[0]])
        logo_plot.set_axis_off()
        f.add_subplot(logo_plot)

        stem_plot = plt.Subplot(f, gs[1], sharex=logo_plot)
        markerline, stemlines, baseline  = stem_plot.stem(X[:fs], [y for y in y_phylop_pixels[:fs]], markerfmt="_", linefmt="-", markerfacecolor=flankingstemcolor, color=greycolor)
        setp(stemlines, 'color', flankingstemcolor)
        setp(markerline, 'markerfacecolor', flankingstemcolor)
        setp(markerline, 'color', flankingstemcolor)
        setp(stemlines, 'linewidth', linewidth)
        setp(markerline, 'markersize', markersize)
        setp(baseline, 'linewidth', linewidth-0.5)
        setp(markerline, 'markeredgewidth', markeredgewidth)
        markerline, stemlines, baseline  = stem_plot.stem(X[fs:-fs], [y for y in y_phylop_pixels[fs:-fs]], markerfmt="g_", linefmt="g-", basefmt="r-")
        setp(stemlines, 'linewidth', linewidth)
        setp(markerline, 'markersize', markersize)
        setp(markerline, 'markeredgewidth', markeredgewidth)
        setp(baseline, 'linewidth', linewidth-0.5)
        markerline, stemlines, baseline  =  stem_plot.stem(X[-fs:], [y for y in y_phylop_pixels[-fs:]], markerfmt="_", linefmt="-", markerfacecolor=flankingstemcolor, color=greycolor)
        setp(stemlines, 'color', flankingstemcolor)
        setp(markerline, 'markerfacecolor', flankingstemcolor)
        setp(stemlines, 'linewidth', linewidth)
        setp(markerline, 'markersize', markersize)
        setp(markerline, 'markeredgewidth', markeredgewidth)
        setp(markerline, 'color', flankingstemcolor)
        setp(baseline, 'linewidth', linewidth-0.5)


        indices_str=[]
        indices1 = np.linspace(-fs,-1, 2)
        for i in indices1:
            indices_str.append('')
        indices2 = np.arange(0, len(X)-2*fs,5)
        for i in indices2:
            indices_str.append('${}$'.format(int(i)+1))

        indices3 = np.linspace(motif_length, motif_length+fs-1, 2)

        for i in indices3:
            indices_str.append('')


        indices12 = np.concatenate((indices1, indices2))
        indices = np.concatenate((indices12, indices3))
        xticks = [X[int(i)+fs] for i in indices]

        max_yticks = 3
        yloc = plt.MaxNLocator(max_yticks)
        stem_plot.yaxis.set_major_locator(yloc)

        #ticks_and_labels = np.linspace(1.02*min(min(y_phylop_pixels), -0.1), 1.02*max(y_phylop_pixels), num = 5, endpoint=True)
        #stem_plot.set_yticks(ticks_and_labels)
        #stem_plot.set_yticklabels(['$%.2f$' %x for x in ticks_and_labels])#(["%0.2f"%(min(y_phylop_pixels)/__scale__), "%0.2f"%(np.mean(y_phylop_pixels)/__scale__), "%0.2f"%(max(y_phylop_pixels)/__scale__)], fontsize=fontsize)
        stem_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Position}$', fontsize=fontsize, fontweight='bold')
        stem_plot.set_xlim([1.2*a,X[-1]+linewidth*1.8])
        stem_plot.set_ylim([min(np.min(y_phylop_pixels), -0.01)-0.03, np.max(y_phylop_pixels,0.01)])
        stem_plot.get_xaxis().tick_bottom()
        stem_plot.get_yaxis().tick_left()
        stem_plot.set_xticks(xticks)
        stem_plot.set_xticklabels(indices_str, fontsize=fontsize)
        stem_plot.spines['top'].set_visible(False)
        stem_plot.spines['right'].set_visible(False)
        stem_plot.yaxis.set_ticks_position('left')
        stem_plot.xaxis.set_ticks_position('bottom')
        stem_plot.spines['left'].set_position('zero')
        #stem_plot.spines['bottom'].set_position(matplotlib.transforms.Bbox(array([[0.125,0.63],[0.25,0.25]])))
        stem_plot.get_yaxis().set_tick_params(direction='out')
        stem_plot.get_xaxis().set_tick_params(direction='out')
        stem_plot.tick_params(axis='y', which='major', pad=tickpad)
        stem_plot.tick_params(axis='x', which='major', pad=tickpad)
        stem_plot.tick_params('both', length=ticklength, width=2, which='major')
        stem_plot.set_ylabel('$\mathrm{PhyloP}\ \mathrm{Score}$', fontsize=fontsize)
        f.add_subplot(stem_plot)

        if sample_gerp_data:
            if annotate:
                gs1 =  gridspec.GridSpec(2, 4, height_ratios=[1,4], width_ratios=[1,1,1,1])
                gerp_header_subplot_gs = gs1[0,1]
                gerp_subplot_gs = gs1[1,1]
                histogram_header_subplot_gs = gs1[0,2]
                histogram_subplot_gs = gs1[1,2]
                ann_header_subplot_gs = gs1[0,3]
                ann_subplot_gs = gs1[1,3]
            else:
                gs1 =  gridspec.GridSpec(2, 3, height_ratios=[1,4], width_ratios=[1,1,1])
                gerp_header_subplot_gs = gs1[0,1]
                gerp_subplot_gs = gs1[1,1]
                histogram_header_subplot_gs = gs1[0,2]
                histogram_subplot_gs = gs1[1,2]
        else:
            if annotate:
                gs1 =  gridspec.GridSpec(2, 3, height_ratios=[1,4], width_ratios=[1,1,1])
                histogram_header_subplot_gs = gs1[0,1]
                histogram_subplot_gs = gs1[1,1]
                ann_header_subplot_gs = gs1[0,2]
                ann_subplot_gs = gs1[1,2]
            else:
                gs1 =  gridspec.GridSpec(2, 2, height_ratios=[1,4], width_ratios=[1,1])
                histogram_header_subplot_gs = gs1[0,1]
                histogram_subplot_gs = gs1[1,1]

        gs1.update(bottom=0.14, right=0.95, left=1-right*0.85, wspace=0.5)


        phlyop_plots_leg = plt.Subplot(f, gs1[0,0], autoscale_on=True)
        pearsonr_pval = str('%.1g'%pr_p[1])
        if 'e' in pearsonr_pval:
            pearsonr_pval += '}'
            pearsonr_pval = pearsonr_pval.replace('e', '*10^{').replace('-0','-')
        score_pval = str('%.1g'%p_deltaphylop)
        if 'e' in score_pval:
            score_pval += '}'
            score_pval = score_pval.replace('e', '*10^{').replace('-0','-')

        textstr = r'\noindent$R_{pearson}=%.2f$($p=%s$)\\~\\$\Delta_{Phylop}=%.2f$($p=%s$)\\~\\' %(pr_p[0], pearsonr_pval, delta_phylop, score_pval)#, reg_phylop_control.rsquared, num_occurrences*reg_phylop_control.params[1])
        txtx = 1-legend_xmultiplier*len(textstr)/100.0
        phlyop_plots_leg.set_frame_on(False)
        phlyop_plots_leg.set_xticks([])
        phlyop_plots_leg.set_yticks([])
        phlyop_plots_leg.text(txtx, txty, textstr, fontsize=legend_fontsize)
        f.add_subplot(phlyop_plots_leg)

        phylop_scatter_plot = plt.Subplot(f, gs1[1,0], autoscale_on=True)
        fit = np.polyfit(motif_scores,motif_sample_phylo_scores,1)
        fit_fn = np.poly1d(fit)

        phylop_scatter_plot.scatter(motif_scores, motif_sample_phylo_scores, color='g', s=[pointsize for i in motif_scores])
        phylop_scatter_plot.plot(motif_scores, y_reg_phylop_sample, 'g', motif_scores, fit_fn(motif_scores), color='g', linewidth=plot_linewidth)
        phylop_scatter_plot.scatter(motif_scores, motif_control_phylo_scores, color=greycolor, s=[pointsize for i in motif_scores])
        phylop_scatter_plot.plot(motif_scores, y_reg_phylop_control, color=greycolor, linewidth=plot_linewidth)

        ticks_and_labels = np.linspace(1.02*min(motif_scores), 1.02*max(motif_scores), num = 5, endpoint=True)
        phylop_scatter_plot.set_xticks(ticks_and_labels)
        ticks_and_labels = ["$%.2f$"%(x/num_occurrences) for x in ticks_and_labels]
        phylop_scatter_plot.set_xticklabels(ticks_and_labels)

        ##max_xticks = 5
        ##xloc = plt.MaxNLocator(max_xticks)
        ##print xloc
        ##phylop_scatter_plot.xaxis.set_major_locator(xloc)
        #ticks_and_labels = np.linspace(1.02*min(min(shifted_sample_phylo_scores), min(shifted_control_phylo_scores)), 1.02*max(max(shifted_sample_phylo_scores),max(shifted_control_phylo_scores)),
                                    #num = 4, endpoint=True)
        #phylop_scatter_plot.set_yticks(ticks_and_labels)
        #phylop_scatter_plot.set_yticklabels(["$%0.2f$"%x for x in ticks_and_labels])
        max_yticks = 4
        yloc = plt.MaxNLocator(max_yticks)
        phylop_scatter_plot.yaxis.set_major_locator(yloc)
        phylop_scatter_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Frequency}$', fontsize=fontsize, fontweight='bold')
        phylop_scatter_plot.get_xaxis().tick_bottom()
        phylop_scatter_plot.get_yaxis().tick_left()
        phylop_scatter_plot.set_ylabel('$\mathrm{PhyloP}\ \mathrm{Score}$', fontsize=fontsize, fontweight='bold')
        phylop_scatter_plot.tick_params(axis='y', which='major', pad=tickpad)
        phylop_scatter_plot.tick_params(axis='x', which='major', pad=tickpad)
        phylop_scatter_plot.get_yaxis().set_tick_params(direction='out')
        phylop_scatter_plot.get_xaxis().set_tick_params(direction='out')
        phylop_scatter_plot.tick_params('both', length=ticklength, width=2, which='major')

        f.add_subplot(phylop_scatter_plot)

        gerp_plots_leg = plt.Subplot(f, gerp_header_subplot_gs, autoscale_on=True)
        gerp_plots_leg.set_frame_on(False)
        gerp_plots_leg.set_xticks([])
        gerp_plots_leg.set_yticks([])
        pearsonr_pval = str('%.1g'%pr_p[1])
        if 'e' in pearsonr_pval:
            pearsonr_pval += '}'
            pearsonr_pval = pearsonr_pval.replace('e', '*10^{').replace('-0','-')

        if sample_gerp_data:
            score_pval = str('%.1g'%p_deltagerp)
            if 'e' in score_pval:
                score_pval += '}'
                score_pval = score_pval.replace('e', '*10^{').replace('-0','-')
            textstr = r'\noindent$R_{pearson}=%.2f$($p=%s$)\\~\\$\Delta_{{Gerp}}=%.2f$($p=%s$)\\~\\'%(pr_g[0], pearsonr_pval, delta_gerp, score_pval)
            txtx = 1-legend_xmultiplier*len(textstr)/100.0
            gerp_plots_leg.text(txtx, txty, textstr, fontsize=legend_fontsize)
            f.add_subplot(gerp_plots_leg)

            gerp_scatter_plot = plt.Subplot(f, gerp_subplot_gs, autoscale_on=True)
            gerp_scatter_plot.scatter(motif_scores, motif_sample_gerp_scores, color='g', s=[pointsize for i in motif_scores])
            gerp_scatter_plot.plot(motif_scores, y_reg_gerp_sample, color='g', linewidth=plot_linewidth)
            gerp_scatter_plot.scatter(motif_scores, motif_control_gerp_scores, color=greycolor, s=[pointsize for i in motif_scores])
            gerp_scatter_plot.plot(motif_scores, y_reg_gerp_control, color=greycolor, linewidth=plot_linewidth)
            ticks_and_labels = np.linspace(1.02*min(motif_scores), 1.02*max(motif_scores), num = 5, endpoint=True)
            gerp_scatter_plot.set_xticks(ticks_and_labels)
            ticks_and_labels = ["$%.2f$"%(x/num_occurrences) for x in ticks_and_labels]
            gerp_scatter_plot.set_xticklabels(ticks_and_labels)

            ##max_xticks = 5
            ##xloc = plt.MaxNLocator(max_xticks)
            ##gerp_scatter_plot.xaxis.set_major_locator(xloc)
            max_yticks = 4
            yloc = plt.MaxNLocator(max_yticks)
            gerp_scatter_plot.yaxis.set_major_locator(yloc)
            gerp_scatter_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Frequency}$', fontsize=fontsize, fontweight='bold')
            gerp_scatter_plot.set_ylabel('$\mathrm{GERP}\ \mathrm{Score}$', fontsize=fontsize, fontweight='bold')
            gerp_scatter_plot.get_xaxis().tick_bottom()
            gerp_scatter_plot.get_yaxis().tick_left()
            gerp_scatter_plot.get_yaxis().set_tick_params(direction='out')
            gerp_scatter_plot.get_xaxis().set_tick_params(direction='out')
            gerp_scatter_plot.tick_params(axis='y', which='major', pad=tickpad)
            gerp_scatter_plot.tick_params(axis='x', which='major', pad=tickpad)
            gerp_scatter_plot.tick_params('both', length=ticklength, width=2, which='major')
            f.add_subplot(gerp_scatter_plot)


        enrichment_plot4 = plt.Subplot(f, histogram_header_subplot_gs, autoscale_on=True)
        enrichment_plot4.set_frame_on(False)
        enrichment_plot4.set_xticks([])
        enrichment_plot4.set_yticks([])
        all_distances = get_motif_distances(peak_file, fimo_file)
        fimo_dir = os.path.dirname(fimo_file)
        motifs_within_100 = filter(lambda x: x<=100 and x>=-100, all_distances)
        motifs_within_100_200 = filter(lambda x: (x<200 and x>100) or (x>-200 and x<-100), all_distances)
        if len(motifs_within_100_200)>0:
            enrichment = len(motifs_within_100)/(len(motifs_within_100_200))#+len(motifs_within_100))
        else:
            enrichment = 1
        enrichment_pval = 0
        number_of_sites = len(motifs_within_100)+len(motifs_within_100_200) #fimo_sites_intersect(fimo_file)
        probability = 200/(ENRICHMENT_SEQ_LENGTH-motif_length)
        enrichment_pval=binom.sf(len(motifs_within_100), number_of_sites, probability)
        enrichment_pval = str('%.1g'%enrichment_pval)
        if 'e' in enrichment_pval:
            enrichment_pval+= '}'
            enrichment_pval= enrichment_pval.replace('e', '*10^{').replace('-0','-')
        textstr = r'\noindent$Enrichment={0:.2f}$\\~\\$(p={1})$'.format(enrichment, enrichment_pval)
        txtx = 0.1*len(textstr)/100.0
        enrichment_plot4.text(txtx, txty, textstr, fontsize=legend_fontsize)
        f.add_subplot(enrichment_plot4)
        enrichment_plot = plt.Subplot(f, histogram_subplot_gs, autoscale_on=True)
        enrichment_plot.hist(all_distances, histogram_nbins, color='white', alpha=0.8, range=[-200,200])
        enrichment_plot.set_xticks([-200,-100,0,100,200])
        max_yticks = 3
        yloc = plt.MaxNLocator(max_yticks)
        enrichment_plot.yaxis.set_major_locator(yloc)
        #enrichment_plot.set_yticks(range(1,6))
        ticks_and_labels = [-200,-100,0,100,200]
        all_distances = np.asarray(all_distances)
        enrichment_plot.set_xticklabels(['${}$'.format(x) for x in ticks_and_labels])
        enrichment_plot.tick_params(axis='y', which='major', pad=tickpad)
        enrichment_plot.tick_params(axis='x', which='major', pad=tickpad)
        enrichment_plot.tick_params('both', length=ticklength, width=2, which='major')
        enrichment_plot.get_xaxis().tick_bottom()
        enrichment_plot.get_yaxis().tick_left()
        enrichment_plot.get_yaxis().set_tick_params(direction='out')
        enrichment_plot.get_xaxis().set_tick_params(direction='out')
        enrichment_plot.axvline(x=-100, linewidth=3, color='red', linestyle='-.')
        enrichment_plot.axvline(x=100, linewidth=3, color='red', linestyle='-.')
        f.add_subplot(enrichment_plot)
        if 'rc' not in ln:
            out_file = os.path.join(fimo_dir,'motif{}Combined_plots.png'.format(motif_number))
            out_file = 'motif{}Combined_plots.png'.format(motif_number)
        else:
            out_file = os.path.join(fimo_dir,'motif{}Combined_plots_rc.png'.format(motif_number))
            out_file = 'motif{}Combined_plots_rc.png'.format(motif_number)

        if annotate:
            filename = r'$'+annotate[0]+'$'
            try:
                a_motif = r'$'+annotate[1]+'$'
            except IndexError:
                a_motif = ''
            try:
                cell_line = r'$'+annotate[2]+'$'
            except IndexError:
                cell_line = ''
            try:
                assay = r'$'+annotate[3]+'$'
            except IndexError:
                assay = ''

            #data = [[r'$Filename$', filename], [r'$Motif$', a_motif], [r'$Cell\ Line$', cell_line], [r'Assay', assay]]
            keys = ['title', 'gene_name', 'dataset', 'assembly']
            data = [[r'$'+key.replace("_", " ").upper()+'$', r'$'+annotate_dict[key]+'$'] for key in keys]
            ann_header = plt.Subplot(f, ann_header_subplot_gs, autoscale_on=True)
            ann_header.set_frame_on(False)
            ann_header.set_xticks([])
            ann_header.set_yticks([])
            f.add_subplot(ann_header)
            textstr = r'$Metadata$'
            txtx = 1.7*len(textstr)/100.0
            ann_header.text(txtx, txty, textstr, fontsize=legend_fontsize)
            ann_plot = plt.Subplot(f, ann_subplot_gs, autoscale_on=True)
            ann_plot.set_xticks([])
            ann_plot.set_yticks([])
            ann_plot.set_frame_on(False)
            table = ann_plot.table(cellText=data,loc='center')
            table.scale(1,2)
            fontproperties=FontProperties(size=legend_fontsize*8)#, family='serif' )
            for key, cell in table.get_celld().items():
                row, col = key
                if row > 0 and col > 0:
                    cell.set_text_props(fontproperties=fontproperties)

            table.set_fontsize(legend_fontsize*8)
            f.add_subplot(ann_plot)

        f.savefig(out_file, figsize=figsize, dpi=dpi)



"""
@click.command()
@click.option('-i', '--meme', help='Meme input file', required=True)
@click.option('-m', '--motif', help='Motif number', default=1, required=True)
@click.option('-ps', '--phylop_sample', help='Sample PhyloP conservation scores', required=True)
@click.option('-pc', '--phylop_control',  help='Control PhyloP conservation scores', required=True)
@click.option('-gs', '--gerp_sample', default=None)
@click.option('-gc', '--gerp_control', default=None)
@click.option('-peak', '--peak_file', help='Path to peaks file')
@click.option('-fimo', '--fimo_file', help='Path to fimo_2_sites output')
@click.option('-f', '--flanking_sites',  default=10)
@click.option('-a', '--annotate', default=None)

def options(meme,motif,flanking_sites,phylop_sample,phylop_control,gerp_sample,gerp_control,peak_file,fimo_file, annotate):
    create_plot(meme,
                motif,
                flanking_sites,
                phylop_sample,
                phylop_control,
                gerp_sample,
                gerp_control,
                peak_file,
                fimo_file,
                annotate)


if __name__ == '__main__':
    options()
"""
