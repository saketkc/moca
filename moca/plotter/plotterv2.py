#!/usr/bin/env python
"""
Generate conservation plots
"""
from __future__ import division
import json
import os
import sys
import matplotlib
matplotlib.use('Agg')
from pylab import setp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import numpy as np

from moca.helpers import read_memefile
from moca.helpers import MocaException

from moca.helpers import get_motif_ic
from moca.helpers import get_motif_bg_freq
from moca.helpers import get_max_occuring_bases
from moca.helpers import create_position_profile

MAGIC_NUM=39.33333333
SCORE_TYPE = 'counts'
# Use 'pssm' or 'counts' to score:
bases = ['A', 'T', 'G', 'C']

## Plot parameters
LINEWIDTH = 3
PHYLOP_BAR_LINEWIDTH = 5
HIST_NBINS = 20


FONTSIZE = 20
LEGEND_FONTSIZE = 26
POINTSIZE = 20
DPI = 300

TICKPAD = 20
TICKLENGTH= 10

GREYNESS = '0.65'
STEM_MARKER_SIZE = 13
STEM_MARKER_EDGEWIDTH = 3.5
STEM_FLANKING_COLOR = 'y'
STEM_LINEWIDTH = 3

##  Text position
TXT_XPOS = 0.02
TXT_YPOS = 0.09

LEGEND_XMULTIPLIER = 1.4
LEGEND_YMULTIPLIER = 1

MAX_YTICKS = 3

def setup_matplotlib():
    """Setup matplotlib
    """
    plt.rc('text', usetex=True)
    plt.rc('font', family='monospace')
    plt.rcParams.update({'axes.titlesize': 'small'})
    plt.rcParams.update({'backend' : 'Agg'})
    plt.rcParams['xtick.labelsize'] = FONTSIZE
    plt.rcParams['ytick.labelsize'] = FONTSIZE
    plt.rcParams['text.latex.preamble'] = [r'\boldmath']

def save_figure(fig, ax, name, form):
    """Save figure to file
    """
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig('{}.{}'.format(name, form), format=form, dpi=DPI, bbox_inches=extent.expanded(1.2, 1.1))
    if form!='png':
        fig.savefig('{}.{}'.format(name, 'png'), format='png', dpi=DPI, bbox_inches=extent.expanded(1.2, 1.1))

def setp_lines(stemlines, markerline, baseline):
    """Set up format for stem plots"""
    setp(stemlines, 'linewidth', STEM_LINEWIDTH)
    setp(markerline, 'markersize', STEM_MARKER_SIZE)
    setp(baseline, 'linewidth', STEM_LINEWIDTH-0.5)
    setp(markerline, 'markeredgewidth', STEM_MARKER_EDGEWIDTH)

def create_stemplot(matplot_dict, X_values, Y_values, motif_length, flank_length=0):
    """Create stem plot for phylop/scores scoresi

    Parameters
    ----------
    matplot_dict: dict like
        A dict like object with the following fields: {'figure': plt.figure, 'gridspec': plt.gridspec, 'shareX': plt.axes}
        where 'gridspec' represents the grid specification  and shareX represents the axis to share X axis with.

    figure: plt.figure
        matplotlib figure

    X_values: np.array
        Array of X values

    Y_values: np.array
        Array of Y-values

    motif_length: int
        Motif length(redundant)

    flank_length: int
        Number of flanking sites
    """
    indices_str=[]
    indices_left = np.linspace(-flank_length, -1, 2)
    indices_str = ['' for x in indices_left]
    indices_center = np.arange(0, len(X_values)-2*flank_length, 5)
    for i in indices_center:
        indices_str.append('${}$'.format(int(i)+1))

    indices_right = np.linspace(motif_length, motif_length+flank_length-1, 2)

    for i in indices_right:
        indices_str.append('')

    indices = np.hstack((indices_left, indices_center, indices_right))
    xticks = [X_values[int(i)+flank_length] for i in indices]
    if len(X_values)!=len(Y_values):
        raise MocaException('Error creating stem plots. (X,Y) dimenstion mismatch:\
                            ({},{})'.format(len(X_values), len(Y_values)))
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    shareX = matplot_dict['shareX']
    stem_plot = plt.Subplot(f, gs, sharex=shareX)
    if flank_length>0:
        X_flank_left = X_values[:flank_length]
        Y_flank_left = Y_values[:flank_length]
        X_center = X_values[flank_length:-flank_length]
        Y_center = Y_values[flank_length:-flank_length]
        X_flank_right = X_values[-flank_length:]
        Y_flank_right = Y_values[-flank_length:]
    else:
        X_center = X_values
        Y_center = Y_values

    markerline, stemlines, baseline  = stem_plot.stem(X_center, Y_center,
                                                      markerfmt="g_",
                                                      linefmt="g-",
                                                      basefmt="r-")
    setp_lines(stemlines, markerline, baseline)
    if flank_length>0:
        markerline, stemlines, baseline  = stem_plot.stem(X_flank_left,
                                                          Y_flank_left,
                                                          markerfmt="_",
                                                          linefmt="-",
                                                          markerfacecolor=STEM_FLANKING_COLOR,
                                                          color=GREYNESS)
        setp(stemlines, 'color', STEM_FLANKING_COLOR)
        setp(markerline, 'markerfacecolor', STEM_FLANKING_COLOR)
        setp_lines(stemlines, markerline, baseline)
        setp(markerline, 'color', STEM_FLANKING_COLOR)
        markerline, stemlines, baseline  =  stem_plot.stem(X_flank_right, Y_flank_right,
                                                           markerfmt="_",
                                                           linefmt="-",
                                                           markerfacecolor=STEM_FLANKING_COLOR,
                                                           color=GREYNESS)
        setp(stemlines, 'color', STEM_FLANKING_COLOR)
        setp(markerline, 'markerfacecolor', STEM_FLANKING_COLOR)
        setp(markerline, 'color', STEM_FLANKING_COLOR)
        setp_lines(stemlines, markerline, baseline)

    yloc = plt.MaxNLocator(MAX_YTICKS)
    stem_plot.yaxis.set_major_locator(yloc)

    stem_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Position}$',
                         fontsize=FONTSIZE,
                         fontweight='bold')
    stem_plot.set_xlim([1.2*MAGIC_NUM, X_values[-1]+LINEWIDTH*1.8])
    stem_plot.set_ylim([min(np.min(Y_values), -0.01)-0.03, np.max(Y_values,0.01)])
    stem_plot.get_xaxis().tick_bottom()
    stem_plot.get_yaxis().tick_left()
    stem_plot.set_xticks(xticks)
    stem_plot.set_xticklabels(indices_str, fontsize=FONTSIZE)
    stem_plot.spines['top'].set_visible(False)
    stem_plot.spines['right'].set_visible(False)
    stem_plot.yaxis.set_ticks_position('left')
    stem_plot.xaxis.set_ticks_position('bottom')
    stem_plot.spines['left'].set_position('zero')
    stem_plot.get_yaxis().set_tick_params(direction='out')
    stem_plot.get_xaxis().set_tick_params(direction='out')
    stem_plot.tick_params(axis='y', which='major', pad=TICKPAD)
    stem_plot.tick_params(axis='x', which='major', pad=TICKPAD)
    stem_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')
    stem_plot.set_ylabel('$\mathrm{PhyloP}\ \mathrm{Score}$', fontsize=FONTSIZE)
    f.add_subplot(stem_plot)

def create_logo_plot(matplot_dict, meme_dir, motif_number, motif_length):
    """Create stem plot for phylop/scores scoresi

    Parameters
    ----------
    matplot_dict: dict like
        A dict like object with the following fields: {'figure': plt.figure, 'gridspec': plt.gridspec, 'shareX': plt.axes}
        where 'gridspec' represents the grid specification  and shareX represents the axis to share X axis with.

    """
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    logo_plot = plt.Subplot(f, gs)
    logo = plt.imread(os.path.join(meme_dir, 'logo{}.png'.format(motif_number)))
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

    logo_plot.imshow(logo, extent=[40+15+z*(MAGIC_NUM+1.9),logo.shape[1]+15+XSCALE_FACTOR*(MAGIC_NUM+1.9),0,logo.shape[0]])
    logo_plot.set_axis_off()
    f.add_subplot(logo_plot)


def create_regression_plot(plot_type='gerp'):
    pass

def _get_logo_path(meme_dir, motif, rc=False):
    return os.path.join(meme_dir, 'logo_rc{}.png'.format(motif) if rc else 'logo{}.png'.format(motif))

def init_figure(meme_dir=None, X_values=None, motif=1, use_gerp=False, annotate=False):
    """Initialize plot based on logo size"""
    logo = plt.imread(_get_logo_path(meme_dir, motif))
    height_px = logo.shape[0] # Should be 212

    if use_gerp:
        if annotate:
            total_px = X_values[-1]+8*height_px+140
            right = (8*height_px+10+140-0.2*height_px)/total_px
        else:
            total_px = X_values[-1]+6*height_px+140
            right = (6*height_px+10+140-0.2*height_px)/total_px
    else:
        if annotate:
            total_px = X_values[-1]+6*height_px+140
            right = (6*height_px+10+140-0.2*height_px)/total_px
        else:
            total_px = X_values[-1]+4*height_px+140
            right = (4*height_px+10+140-0.2*height_px)/total_px

    figsize=(total_px/100,(2*height_px)/100+0.6)

    gs =  gridspec.GridSpec(2, 1)#, width_ratios=[1, right], height_ratios=[1,1])
    gs.update(top=1.0, bottom=0.14, left=0.08, right=1-right)#, right=0.8)#, left=0.06)#, right=right, wspace=0.025, hspace=0.03, wd)
    f = plt.figure(figsize=figsize, dpi=DPI, facecolor='w', edgecolor='k')

    return {'figure': f,
            'gs': gs,
            'figsize': figsize,
            'right_margin': right,
            'total_px': total_px}

def enrichment_plot(matplot_dict, motif_length, peak_file, fimo_file):
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    shareX = matplot_dict['shareX']
    enrichment_plot = plt.Subplot(f, gs, autoscale_on=True)
    enrichment_plot.set_frame_on(False)
    enrichment_plot.set_xticks([])
    enrichment_plot.set_yticks([])
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
    enrichment_plot.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
    f.add_subplot(enrichment_plot)
    enrichment_plot = plt.Subplot(f, histogram_subplot_gs, autoscale_on=True)
    enrichment_plot.hist(all_distances, histogram_nbins, color='white', alpha=0.8, range=[-200,200])
    enrichment_plot.set_xticks([-200,-100,0,100,200])
    max_yticks = 3
    yloc = plt.MaxNLocator(max_yticks)
    enrichment_plot.yaxis.set_major_locator(yloc)
    ticks_and_labels = [-200,-100,0,100,200]
    all_distances = np.asarray(all_distances)
    enrichment_plot.set_xticklabels(['${}$'.format(x) for x in ticks_and_labels])
    enrichment_plot.tick_params(axis='y', which='major', pad=TICKPAD)
    enrichment_plot.tick_params(axis='x', which='major', pad=TICKPAD)
    enrichment_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')
    enrichment_plot.get_xaxis().tick_bottom()
    enrichment_plot.get_yaxis().tick_left()
    enrichment_plot.get_yaxis().set_tick_params(direction='out')
    enrichment_plot.get_xaxis().set_tick_params(direction='out')
    enrichment_plot.axvline(x=-100, linewidth=3, color='red', linestyle='-.')
    enrichment_plot.axvline(x=100, linewidth=3, color='red', linestyle='-.')
    f.add_subplot(enrichment_plot)

def create_annnotation_plot(matplot_dict, annotate):
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    shareX = matplot_dict['shareX']
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
    keys = ['title', 'gene_name', 'dataset', 'assembly']
    data = [[r'$'+key.replace("_", " ").upper()+'$', r'$'+annotate_dict[key]+'$'] for key in keys]
    ann_header = plt.Subplot(f, ann_header_subplot_gs, autoscale_on=True)
    ann_header.set_frame_on(False)
    ann_header.set_xticks([])
    ann_header.set_yticks([])
    f.add_subplot(ann_header)
    textstr = r'$Metadata$'
    txtx = 1.7*len(textstr)/100.0
    ann_header.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
    ann_plot = plt.Subplot(f, ann_subplot_gs, autoscale_on=True)
    ann_plot.set_xticks([])
    ann_plot.set_yticks([])
    ann_plot.set_frame_on(False)
    table = ann_plot.table(cellText=data,loc='center')
    table.scale(1,2)
    fontproperties=FontProperties(size=LEGEND_FONTSIZE*8)#, family='serif' )
    for key, cell in table.get_celld().items():
        row, col = key
        if row > 0 and col > 0:
            cell.set_text_props(fontproperties=fontproperties)

    table.set_fontsize(LEGEND_FONTSIZE*8)
    f.add_subplot(ann_plot)

def create_phylop_leg_plot(matplot_dict):
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    shareX = matplot_dict['shareX']

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
    txtx = 1-LEGEND_XMULTIPLIER*len(textstr)/100.0
    phlyop_plots_leg.set_frame_on(False)
    phlyop_plots_leg.set_xticks([])
    phlyop_plots_leg.set_yticks([])
    phlyop_plots_leg.text(txtx, txty, textstr, fontsize=LEGEND_FONTSIZE)
    f.add_subplot(phlyop_plots_leg)

def create_phylop_scatter(matplot_dict, motif_freq):

    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    phylop_scatter_plot = plt.Subplot(f, gs1[1,0], autoscale_on=True)
    fit = np.polyfit(motif_freq,motif_sample_phylo_scores,1)
    fit_fn = np.poly1d(fit)

    phylop_scatter_plot.scatter(motif_freq, motif_sample_phylo_scores, color='g', s=[pointsize for i in motif_freq])
    phylop_scatter_plot.plot(motif_freq, y_reg_phylop_sample, 'g', motif_freq, fit_fn(motif_freq), color='g', linewidth=plot_linewidth)
    phylop_scatter_plot.scatter(motif_freq, motif_control_phylo_scores, color=greycolor, s=[pointsize for i in motif_freq])
    phylop_scatter_plot.plot(motif_freq, y_reg_phylop_control, color=greycolor, linewidth=plot_linewidth)

    ticks_and_labels = np.linspace(1.02*min(motif_freq), 1.02*max(motif_freq), num = 5, endpoint=True)
    phylop_scatter_plot.set_xticks(ticks_and_labels)
    ticks_and_labels = ["$%.2f$"%(x/num_occurrences) for x in ticks_and_labels]
    phylop_scatter_plot.set_xticklabels(ticks_and_labels)

    max_yticks = 4
    yloc = plt.MaxNLocator(max_yticks)
    phylop_scatter_plot.yaxis.set_major_locator(yloc)
    phylop_scatter_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Frequency}$', fontsize=fontsize, fontweight='bold')
    phylop_scatter_plot.get_xaxis().tick_bottom()
    phylop_scatter_plot.get_yaxis().tick_left()
    phylop_scatter_plot.set_ylabel('$\mathrm{PhyloP}\ \mathrm{Score}$', fontsize=fontsize, fontweight='bold')
    phylop_scatter_plot.tick_params(axis='y', which='major', pad=TICKPAD)
    phylop_scatter_plot.tick_params(axis='x', which='major', pad=TICKPAD)
    phylop_scatter_plot.get_yaxis().set_tick_params(direction='out')
    phylop_scatter_plot.get_xaxis().set_tick_params(direction='out')
    phylop_scatter_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')

    f.add_subplot(phylop_scatter_plot)

def create_gerp_plot(matplot_dict):
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    shareX = matplot_dict['shareX']
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
        txtx = 1-LEGEND_XMULTIPLIER*len(textstr)/100.0
        gerp_plots_leg.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
        f.add_subplot(gerp_plots_leg)

        gerp_scatter_plot = plt.Subplot(f, gerp_subplot_gs, autoscale_on=True)
        gerp_scatter_plot.scatter(motif_freq, motif_sample_gerp_scores, color='g', s=[pointsize for i in motif_freq])
        gerp_scatter_plot.plot(motif_freq, y_reg_gerp_sample, color='g', linewidth=plot_linewidth)
        gerp_scatter_plot.scatter(motif_freq, motif_control_gerp_scores, color=greycolor, s=[pointsize for i in motif_freq])
        gerp_scatter_plot.plot(motif_freq, y_reg_gerp_control, color=greycolor, linewidth=plot_linewidth)
        ticks_and_labels = np.linspace(1.02*min(motif_freq), 1.02*max(motif_freq), num = 5, endpoint=True)
        gerp_scatter_plot.set_xticks(ticks_and_labels)
        ticks_and_labels = ["$%.2f$"%(x/num_occurrences) for x in ticks_and_labels]
        gerp_scatter_plot.set_xticklabels(ticks_and_labels)

        max_yticks = 4
        yloc = plt.MaxNLocator(max_yticks)
        gerp_scatter_plot.yaxis.set_major_locator(yloc)
        gerp_scatter_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Frequency}$', fontsize=fontsize, fontweight='bold')
        gerp_scatter_plot.set_ylabel('$\mathrm{GERP}\ \mathrm{Score}$', fontsize=fontsize, fontweight='bold')
        gerp_scatter_plot.get_xaxis().tick_bottom()
        gerp_scatter_plot.get_yaxis().tick_left()
        gerp_scatter_plot.get_yaxis().set_tick_params(direction='out')
        gerp_scatter_plot.get_xaxis().set_tick_params(direction='out')
        gerp_scatter_plot.tick_params(axis='y', which='major', pad=TICKPAD)
        gerp_scatter_plot.tick_params(axis='x', which='major', pad=TICKPAD)
        gerp_scatter_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')
        f.add_subplot(gerp_scatter_plot)

def create_plot(meme_file,
                peak_file,
                fimo_file,
                sample_phylop_file,
                control_phylop_file,
                sample_gerp_file=None,
                control_gerp_file=None,
                annotate=None,
                motif_number=1,
                flank_length=5):
    meme_record = read_memefile(meme_file)
    record = meme_record['motif_records'][motif_number-1]
    motif_lenth = record.length
    num_occurrences = getattr(record, 'num_occurrences', 'Unknown')
    meme_dir = os.path.abspath(os.path.dirname(meme_file))
    fimo_dir = os.path.abspath(os.path.dirname(fimo_file))

    annotate_dict = None

    if annotate == "" or annotate == ' ':
        annotate = None
    elif annotate:
        with open(annotate) as f:
            annotate_dict = json.load(f)

    if sample_gerp_file:
        use_gerp = False
    else:
        use_gerp = True

    #profile = create_position_profile(record, SCORE_TYPE)
    max_occur = get_max_occuring_bases(record, max_count=1, count_type=SCORE_TYPE)
    ## motif_freq is tn array of scores of the max  occuring base at each position of the motif
    motif_freq = []
    for position in max_occur:
        motif_freq.append(position[0][1])

    motif_freq = np.asarray(motif_freq)
    sample_phylo_scores = np.loadtxt(sample_phylop_file)

    motif = record
    motif_length = motif.length
    meme_dir = os.path.abspath(os.path.dirname(meme_file))
    X = [40+15] ## this is by trial and error, the position for the first base logo
    ## Generate all other X coordinates
    for j in range(1,len(motif)+2*flank_length):
        X.append( X[j-1]+MAGIC_NUM+1.9 )

    motif_bits = get_motif_ic(meme_file, motif_number)

    ##FIXME This is a big dirty hacl to get thegenerate plots for the Reverse complement logo too
    logo_name =['logo{}.png'.format(motif_number), 'logo_rc{}.png'.format(motif_number)]
    for ln in logo_name:
        if 'rc'in ln:
            sample_phylo_scores = sample_phylo_scores[::-1]
        matplot_dict =init_figure(meme_dir=meme_dir, X_values=X, motif=motif_number,
                                    use_gerp=use_gerp, annotate=annotate)
        f = matplot_dict['figure']
        gs = matplot_dict['gs']
        figsize = matplot_dict['figsize']
        right_margin = matplot_dict['right_margin']
        total_px= matplot_dict['total_px']

        create_logo_plot({'figure':f, 'gridspec': gs[0]}, meme_dir, motif_number, motif_length)


        if use_gerp:
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

        gs1.update(bottom=0.14, right=0.95, left=1-right_margin*0.85, wspace=0.5)

        #create_phylop_plot(f,gs1[0,0])

        #create_phylop_scatter(f, gs[1,0])



        if 'rc' not in ln:
            out_file = os.path.join(fimo_dir,'motif{}Combined_plots.png'.format(motif_number))
            out_file = 'motif{}Combined_plots.png'.format(motif_number)
        else:
            out_file = os.path.join(fimo_dir,'motif{}Combined_plots_rc.png'.format(motif_number))
            out_file = 'motif{}Combined_plots_rc.png'.format(motif_number)

        if annotate:
            cretate_annnotation_plot()

        f.savefig(out_file, figsize=figsize, dpi=DPI)

if __name__ == '__main__':
    meme_file, peak_file, fimo_file, sample_phylop_file, control_phylop_file, sample_gerp_file,  control_gerp_file= sys.argv[1:]
    create_plot(meme_file,
                peak_file,
                fimo_file,
                sample_phylop_file,
                control_phylop_file,
                sample_gerp_file=None,
                control_gerp_file=None,
                annotate=None,
                motif_number=1,
                flank_length=5)
