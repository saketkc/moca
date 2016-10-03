#!/usr/bin/anv python
"""
Generate conservation plots
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import six
from builtins import zip
from builtins import str
from builtins import range
import click
import json
import math
import matplotlib
matplotlib.use('Agg')
import os
import numpy as np
import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from pylab import setp
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
click.disable_unicode_literals_warning = True

from moca.helpers import get_max_occuring_bases
from moca.helpers import get_total_sequences
from moca.helpers import MocaException
from moca.helpers import read_centrimo_stats
from moca.helpers import read_centrimo_txt
from moca.helpers import read_memefile
from moca.helpers import safe_makedir
from moca.helpers.seqstats import format_pvalue
from moca.helpers.seqstats import get_flanking_scores
from moca.helpers.seqstats import get_pearson_corr
from moca.helpers.seqstats import perform_OLS
from moca.helpers.seqstats import perform_t_test
from moca.helpers.seqstats import remove_flanking_scores

OFFSET = 39.33333333
COUNT_TYPE = 'counts'
# Use 'pssm' or 'counts' to score:
bases = ['A', 'T', 'G', 'C']

## Plot parameters
LINEWIDTH = 3
PHYLOP_BAR_LINEWIDTH = 5
HIST_NBINS = 20

FONTSIZE = 20
LEGEND_FONTSIZE = 26
POINTSIZE = 24
DPI = 300

TICKPAD = 20
TICKLENGTH= 10

GREYNESS = '0.55'
STEM_MARKER_SIZE = 13
STEM_MARKER_EDGEWIDTH = 3.5
STEM_FLANKING_COLOR = 'y'
STEM_LINEWIDTH = 3

##  Text position
TXT_XPOS = 0.02
TXT_YPOS = 0.09

LEGEND_XMULTIPLIER = 1.8
LEGEND_YMULTIPLIER = 1

MAX_YTICKS = 4

def setup_matplotlib():
    """Setup matplotlib
    """
    plt.rc('text', usetex=True)
    plt.rc('font', family='monospace', weight='bold')
    plt.rcParams.update({'axes.titlesize': 'small'})
    plt.rcParams.update({'backend' : 'Agg'})
    plt.rcParams['xtick.labelsize'] = FONTSIZE
    plt.rcParams['ytick.labelsize'] = FONTSIZE
    plt.rcParams['text.latex.preamble'] = [r'\boldmath']
    plt.rcParams['legend.loc'] = 'best'

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

def create_stemplot(matplot_dict, X_values, Y_values,
                    motif_length, legend_title, flank_length=0):
    """Create stem plot for phylop/scores scoresi

    Parameters
    ----------
    matplot_dict: dict like
        A dict like object with the following fields:
            {'figure': plt.figure, 'gridspec': plt.gridspec, 'shareX': plt.axes}
        where 'gridspec' represents the grid specification
        and shareX represents the axis to share X axis with.

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
    if not len(Y_values):
        return
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
    #if len(X_values)!=len(Y_values):
    #    raise MocaException('Error creating stem plots. (X,Y) dimenstion mismatch:\
    #                        ({},{})'.format(len(X_values), len(Y_values)))
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
    stem_plot.set_xticks(xticks)
    stem_plot.set_xticklabels(indices_str, fontsize=FONTSIZE)
    seaborn.despine(ax=stem_plot, offset=10, trim=True)
    stem_plot.set_ylabel('$\mathrm{%s}\ \mathrm{Score}$'%(legend_title), fontsize=FONTSIZE)
    f.add_subplot(stem_plot)
    return X_flank_left, X_center, X_flank_right

def create_logo_plot(matplot_dict, meme_dir, logo_path, motif_length):
    """Create stem plot for phylop/scores scoresi

    Parameters
    ----------
    matplot_dict: dict like
        A dict like object with the following fields:
            {'figure': plt.figure, 'gridspec': plt.gridspec, 'shareX': plt.axes}
        where 'gridspec' represents the grid specification
        and shareX represents the axis to share X axis with.
    """
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    logo_plot = plt.Subplot(f, gs)
    logo = plt.imread(os.path.join(os.path.abspath(meme_dir), logo_path))
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

    logo_plot.imshow(logo, extent=[40+15+z*(OFFSET+1.9),logo.shape[1]+15+XSCALE_FACTOR*(OFFSET+1.9),0,logo.shape[0]])
    logo_plot.set_axis_off()
    f.add_subplot(logo_plot)
    return logo_plot

def create_bar_plot(logo_plot,  X_right, height_px,
                    total_sequences, all_meme_occurrences,
                    motif_number, motif_evalue):
    start_point = len(all_meme_occurrences)
    height_scale = 2.25
    bottom_scale = 3
    heights = np.array([height_px/height_scale*all_meme_occurrences[i]/total_sequences for i in range(0, start_point)])
    bottoms = np.array([height_px/bottom_scale for i in range(0, start_point)])
    barlist = logo_plot.bar(np.array(X_right[-start_point:]),
                            heights,
                            width=24,
                            bottom=bottoms,
                            fill=False,
                            edgecolor='black')
    occurrence = all_meme_occurrences[motif_number-1]/total_sequences*100.0
    textstr = r'\noindent$\mathrm{E-Value}=%s$\\~\\$\mathrm{Enrichment}=%s$'%(format_pvalue(motif_evalue), '{}\%'.format(occurrence))
    index = int(math.floor(-start_point/2))
    total_heights = heights+bottoms
    max_height = np.max(total_heights)
    logo_plot.text(X_right[0], bottoms[0]/2.5, textstr, fontsize=12)
    barlist[motif_number-1].set_color('red')
    barlist[motif_number-1].set_hatch('/')

def _get_logo_path(meme_dir, motif, rc=False):
    return os.path.join(meme_dir, 'logo_rc{}.png'.format(motif) if rc else 'logo{}.png'.format(motif))

def init_figure(meme_dir=None, X_values=None,
                motif=1, subplot_ncols=None,
                annotate=False):
    """Initialize plot based on logo size"""
    logo = plt.imread(_get_logo_path(meme_dir, motif))
    height_px = logo.shape[0] # Should be 212
    #subplot_ncols+=1
    total_px = X_values[-1]+2*subplot_ncols*height_px+140
    right = (2*subplot_ncols*height_px+10+140-0.2*height_px)/total_px

    figsize=(total_px/100,(2*height_px)/100+0.6)

    gs =  gridspec.GridSpec(2, 1)
    gs.update(top=1.0, bottom=0.14,
              left=0.08, right=1-right)
    f = plt.figure(figsize=figsize, dpi=DPI,
                   facecolor='w', edgecolor='k')

    return {'figure': f,
            'gs': gs,
            'figsize': figsize,
            'right_margin': right,
            'height_px': height_px,
            'total_px': total_px}

def create_enrichment_plot(matplot_dict, motif_number, centrimo_txt, centrimo_stats):
    f = matplot_dict['figure']
    gs_h = matplot_dict['gridspec_header']
    gs_b = matplot_dict['gridspec_body']

    all_stats = read_centrimo_stats(centrimo_stats)
    motif_stats = all_stats['MEME']['MOTIF_{}'.format(motif_number)]
    X_values = np.array(motif_stats['pos'])
    Y_values = np.array(motif_stats['count'])
    normalized_Y = Y_values/np.sum(Y_values)

    density = gaussian_kde(Y_values)
    x_smooth = np.linspace(X_values.min(), X_values.max(), 100)
    #y_smooth = spline(X_values, normalized_Y, x_smooth)
    smoother = UnivariateSpline(X_values, normalized_Y, s=0.0005)
    density.covariance_factor = lambda : .0005
    density._compute_covariance()


    y_smooth = smoother(x_smooth)

    centrimo_dict = read_centrimo_txt(centrimo_txt)
    try:
        centrimo_motif_dict = centrimo_dict[motif_number-1]
    except IndexError:
        #TODO This should be logged somewhere
        centrimo_motif_dict = centrimo_dict[0]

    enrichment_pval = float(centrimo_motif_dict['adj_p-value'])
    enrichment = float(centrimo_motif_dict['sites_in_bin'])/float(centrimo_motif_dict['total_sites'])

    enrichment_plot = plt.Subplot(f, gs_h, autoscale_on=True)
    enrichment_plot.set_frame_on(False)
    enrichment_plot.set_xticks([])
    enrichment_plot.set_yticks([])

    enrichment_pval = str('%.1g'%enrichment_pval)
    if 'e' in enrichment_pval:
        enrichment_pval+= '}'
        enrichment_pval= enrichment_pval.replace('e', '*10^{').replace('-0','-')
    text = '$\mathrm{Center-enrichment}$'
    textstr = r'\noindent {}=${:.2f}$\\~\\$(p={})$'.format(text, enrichment, enrichment_pval)
    txtx = -0.05*len(textstr)/100.0
    enrichment_plot.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE-2)
    f.add_subplot(enrichment_plot)
    enrichment_plot = plt.Subplot(f, gs_b, autoscale_on=True)

    #enrichment_plot.plot(X_values, normalized_Y, linewidth=LINEWIDTH)
    #enrichment_plot.plot(x_smooth, density(x_smooth), linewidth=LINEWIDTH)
    enrichment_plot.plot(x_smooth, y_smooth, linewidth=LINEWIDTH)
    enrichment_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')
    enrichment_plot.set_xlabel('$\mathrm{Distance}\ \mathrm{from} \ \mathrm{peak}$', fontsize=FONTSIZE, fontweight='bold')
    enrichment_plot.set_ylabel('$\mathrm{Density}$', fontsize=FONTSIZE, fontweight='bold')
    enrichment_plot.axvline(x=X_values.min()/2, linewidth=3, color='green', linestyle='-.')
    enrichment_plot.axvline(x=X_values.max()/2, linewidth=3, color='green', linestyle='-.')
    enrichment_plot.axvline(x=0, linewidth=3, color='green', linestyle='-.')
    seaborn.despine(ax=enrichment_plot, offset=10, trim=True)

    enrichment_plot.axvline(x=X_values.min(), linewidth=3, color='red', linestyle='-.')
    enrichment_plot.axvline(x=X_values.max(), linewidth=3, color='red', linestyle='-.')
    MAX_YTICKS = 3
    yloc = plt.MaxNLocator(MAX_YTICKS)
    enrichment_plot.yaxis.set_major_locator(yloc)
    f.add_subplot(enrichment_plot)

def create_annnotation_plot(matplot_dict, json_annotation):
    f = matplot_dict['figure']
    gs_h = matplot_dict['gridspec_header']
    gs_b = matplot_dict['gridspec_body']

    annotate_dict = None
    with open(json_annotation) as f:
        annotate_dict = json.load(f)
    keys = ['title', 'gene_name', 'dataset', 'assembly', 'filename']
    data = [[r'$'+key.replace("_", " ").upper()+'$', r'$'+annotate_dict[key]+'$'] for key in keys]
    ann_header = plt.Subplot(f, gs_h, autoscale_on=True)
    ann_header.set_frame_on(False)
    ann_header.set_xticks([])
    ann_header.set_yticks([])
    f.add_subplot(ann_header)
    textstr = r'$Metadata$'
    txtx = 1.7*len(textstr)/100.0
    ann_header.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
    ann_plot = plt.Subplot(f, gs_b, autoscale_on=True)
    ann_plot.set_xticks([])
    ann_plot.set_yticks([])
    ann_plot.set_frame_on(False)
    table = ann_plot.table(cellText=data,loc='center')
    table.scale(1,2)
    fontproperties=FontProperties(size=LEGEND_FONTSIZE*8)#, family='serif' )
    for key, cell in list(table.get_celld().items()):
        row, col = key
        if row > 0 and col > 0:
            cell.set_text_props(fontproperties=fontproperties)

    table.set_fontsize(LEGEND_FONTSIZE*8)
    f.add_subplot(ann_plot)

def create_ols_legend_plot(matplot_dict, motif_freq,
                           sample_scores, control_scores,
                           flank_length, legend_title):
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']

    phlyop_plots_legend = plt.Subplot(f, gs, autoscale_on=True)
    corr_result = get_pearson_corr(motif_freq,
                                   remove_flanking_scores(sample_scores, flank_length))
    corr_pval = corr_result[1]
    corr_r2 = corr_result[0]
    ttest_result = perform_t_test(remove_flanking_scores(sample_scores, flank_length),
                                  get_flanking_scores(sample_scores, flank_length))
    p_deltaphylop = ttest_result['one_sided_pval']
    delta_phylop = ttest_result['delta']
    #T_deltaphylop = ttest_result['T']
    pearsonr_pval = str('%.1g'%corr_pval)
    if 'e' in pearsonr_pval:
        pearsonr_pval += '}'
        pearsonr_pval = pearsonr_pval.replace('e', '*10^{').replace('-0','-')
    score_pval = str('%.1g'%p_deltaphylop)
    if 'e' in score_pval:
        score_pval += '}'
        score_pval = score_pval.replace('e', '*10^{').replace('-0','-')

    textstr = r'$r^2_{pearson}=%.2f(p=%s)$' '\n' r'$\Delta_{%s}=%.2f(p=%s)$' %(corr_r2, pearsonr_pval, legend_title, delta_phylop, score_pval)
    txtx = 1-LEGEND_XMULTIPLIER*len(textstr)/100.0
    phlyop_plots_legend.set_frame_on(False)
    phlyop_plots_legend.set_xticks([])
    phlyop_plots_legend.set_yticks([])
    phlyop_plots_legend.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
    f.add_subplot(phlyop_plots_legend)


def create_scatter_plot(matplot_dict, motif_freq,
                        sample_scores, control_scores,
                        flank_length, num_occurrences, y_label):

    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    scatter_plot = plt.Subplot(f, gs, autoscale_on=True)
    control_scores = remove_flanking_scores(control_scores, flank_length)
    sample_scores = remove_flanking_scores(sample_scores, flank_length)

    fit = np.polyfit(motif_freq, sample_scores, 1)
    fit_fn = np.poly1d(fit)

    control_ols = perform_OLS(control_scores, motif_freq)
    sample_ols = perform_OLS(sample_scores, motif_freq)

    sample_regression_line = sample_ols['regression_line']
    control_regression_line = control_ols['regression_line']

    s1 = scatter_plot.scatter(motif_freq, sample_scores, color='g',
                              s=[POINTSIZE for i in motif_freq],
                              marker='^', label=r'$\mathrm{Sample}$')
    scatter_plot.plot(motif_freq, sample_regression_line, 'g',
                      motif_freq, fit_fn(motif_freq),
                      color='g', linewidth=LINEWIDTH)
    s2 = scatter_plot.scatter(motif_freq, control_scores,
                              color=GREYNESS, s=[POINTSIZE for i in motif_freq],
                              marker='o', label=r'$\mathrm{Control}$')
    scatter_plot.plot(motif_freq, control_regression_line,
                      color=GREYNESS, linewidth=LINEWIDTH)
    leg = scatter_plot.legend(fontsize=14)
    leg.draw_frame(True)
    #leg.get_frame().set_edgecolor('b')
    leg.get_frame().set_linewidth(2.0)

    ticks_and_labels = np.linspace(1.02*min(motif_freq), 1.02*max(motif_freq),
                                   num = 3, endpoint=True)
    scatter_plot.set_xticks(ticks_and_labels)

    ticks_and_labels = ["$%.2f$"%(x/(1.02*num_occurrences)) for x in ticks_and_labels]
    scatter_plot.set_xticklabels(ticks_and_labels)#, rotation=45)

    yloc = plt.MaxNLocator(MAX_YTICKS)
    scatter_plot.yaxis.set_major_locator(yloc)
    scatter_plot.set_xlabel(r'$\mathrm{Most}\ \mathrm{frequent} \ \mathrm{base}\ \mathrm{frequency}$',
                            fontsize=FONTSIZE, fontweight='bold')
    scatter_plot.get_xaxis().tick_bottom()
    scatter_plot.get_yaxis().tick_left()
    scatter_plot.set_ylabel('$\mathrm{%s}\ \mathrm{Score}$'%(y_label), fontsize=FONTSIZE, fontweight='bold')
    scatter_plot.tick_params(axis='y', which='major', pad=TICKPAD)
    scatter_plot.tick_params(axis='x', which='major', pad=TICKPAD)
    scatter_plot.get_yaxis().set_tick_params(direction='out')
    scatter_plot.get_xaxis().set_tick_params(direction='out')
    scatter_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')

    f.add_subplot(scatter_plot)


def create_plot(meme_file,
                plot_title,
                output_dir=None,
                centrimo_dir=None,
                motif_number=1,
                flank_length=5,
                sample_score_files=[],
                control_score_files=[],
                reg_plot_titles=[],
                annotate=None,
                save=True):
    """Create plot
    Parameters
    ----------
    meme_file: string
        Path to meme.txt
    peak_file: string
        Path to summit file
    centrimo_dir: string
        Path to centrimo's output directory
    motif_number: int
        1-based number of motif in the motif file
    sample_score_files: list
        Path to conservation scores files for sample
    control_score_files: list
        Path to conservation score files for control
    legend_titles: list
        List of legend titles
    """
    meme_record = read_memefile(meme_file)
    total_sequences = get_total_sequences(meme_file)
    record = meme_record['motif_records'][motif_number-1]
    num_occurrences = getattr(record, 'num_occurrences', 'Unknown')
    all_meme_occurrences = []
    for motif_record in meme_record['motif_records']:
        all_meme_occurrences.append(getattr(motif_record, 'num_occurrences', 'Unknown'))

    meme_dir = os.path.abspath(os.path.dirname(meme_file))
    if not output_dir:
        output_dir = os.path.join(os.path.join(meme_dir, '..'), 'moca_plots')
    safe_makedir(output_dir)

    subplot_ncols = 1

    if len(sample_score_files) == 0:
        raise MocaException('Found no sample score files')
    elif len(control_score_files) == 0:
        raise MocaException('Found no control score filees')
    elif len(sample_score_files)!=len(control_score_files):
        raise MocaException('Found unequal size of sample and control score files')

    if annotate == "" or annotate == ' ':
        annotate = None
        subplot_ncols +=1

    max_occur = get_max_occuring_bases(record, max_count=1, count_type=COUNT_TYPE)
    motif_freq = []
    for position in max_occur:
        motif_freq.append(position[0][1])

    motif_freq = np.asarray(motif_freq)
    sample_conservation_scores = []
    control_conservation_scores = []
    for i in range(0, len(sample_score_files)):
        sample_conservation_scores.append(np.loadtxt(sample_score_files[i]))
    for i in range(0, len(control_score_files)):
        control_conservation_scores.append(np.loadtxt(control_score_files[i]))

    motif = record
    motif_length = motif.length
    motif_evalue = motif.evalue
    meme_dir = os.path.abspath(os.path.dirname(meme_file))
    X_values = [40+15] ## this is by trial and error, the position for the first base logo
    ## Generate all other X coordinates
    for j in range(1,len(motif)+2*flank_length):
        X_values.append( X_values[j-1]+OFFSET+1.9 )

    if centrimo_dir:
        subplot_ncols +=1
        centrimo_dir = os.path.abspath(centrimo_dir)
        centrimo_txt = os.path.join(centrimo_dir, 'centrimo.txt')
        centrimo_stats = os.path.join(centrimo_dir, 'site_counts.txt')

    plot_title += r' \# {}'.format(motif_number)
    ##FIXME This is a big dirty hacl to get thegenerate plots for the Reverse complement logo too
    logo_name =['logo{}.png'.format(motif_number), 'logo_rc{}.png'.format(motif_number)]
    figures = []
    for sample_score, control_score, subplot_legend_title in zip(sample_conservation_scores,
                                                  control_conservation_scores,
                                                  reg_plot_titles):
        for logo_filename in logo_name:
            setup_matplotlib()
            if 'rc'in logo_filename:
                sample_score = sample_score[::-1]
            matplot_dict = init_figure(meme_dir=meme_dir, X_values=X_values,
                                    motif=motif_number,
                                    subplot_ncols=subplot_ncols, annotate=annotate)
            f = matplot_dict['figure']
            gs = matplot_dict['gs']
            figsize = matplot_dict['figsize']
            right_margin = matplot_dict['right_margin']
            #total_px= matplot_dict['total_px']

            title = r'\textbf{' + '\\underline{'+'{}'.format(plot_title)+'}}'
            f.suptitle(title, fontsize=LEGEND_FONTSIZE)
            logo_plot = create_logo_plot({'figure':f, 'gridspec': gs[0]}, meme_dir, logo_filename, motif_length)

            subgrid = gridspec.GridSpec(2, subplot_ncols, height_ratios=[1,2], width_ratios=[1]*subplot_ncols)
            subgrid.update(bottom=0.14, right=0.9, left=1-right_margin*0.85, wspace=0.58)
            X_left, X_center, X_right = create_stemplot({'figure': f,
                                                        'gridspec': gs[1],
                                                        'shareX': logo_plot},
                                                        X_values,
                                                        sample_score,
                                                        motif_length,
                                                        flank_length=flank_length,
                                                        legend_title=subplot_legend_title)

            create_bar_plot(logo_plot,  X_right, matplot_dict['height_px'],
                            total_sequences, all_meme_occurrences, motif_number, motif_evalue)
            create_ols_legend_plot({'figure':f, 'gridspec': subgrid[0,0]},  motif_freq,
                                sample_score, control_score,
                                flank_length, legend_title=subplot_legend_title)
            create_scatter_plot({'figure':f, 'gridspec': subgrid[1,0]}, motif_freq,
                                sample_score, control_score,
                                flank_length, num_occurrences, y_label=subplot_legend_title)

            if centrimo_dir:
                create_enrichment_plot({'figure': f,
                                        'gridspec_header': subgrid[0,1],
                                        'gridspec_body': subgrid[1,1]},
                                        motif_number,
                                        centrimo_txt,
                                        centrimo_stats)

            if 'rc' not in logo_filename:
                out_file = os.path.join(output_dir,'moca_{}_{}.png'.format(subplot_legend_title, motif_number))
            else:
                out_file = os.path.join(output_dir,'moca_{}_{}_rc.png'.format(subplot_legend_title, motif_number))

            if annotate:
                create_annnotation_plot({'figure': f,
                                        'gridspec_header': subgrid[0,-1],
                                        'gridspec_body': subgrid[1,-1]},
                                        annotate)

            if save:
                f.savefig(out_file, figsize=figsize, dpi=DPI)
            figures.append(f)
            plt.close('all')
    return figures
if __name__ == '__main__':
    create_plot()
