#!/usr/bin/anv python
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
import statsmodels.api as sm

from moca.helpers import read_memefile
from moca.helpers import get_max_occuring_bases
from moca.helpers import read_centrimo_txt
from moca.helpers import read_centrimo_stats

from moca.plotter import perform_t_test
from moca.plotter import get_pearson_corr

MAGIC_NUM=39.33333333
COUNT_TYPE = 'counts'
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

LEGEND_XMULTIPLIER = 1.8
LEGEND_YMULTIPLIER = 1

MAX_YTICKS = 3

def remove_flanking_scores(all_scores, flank_length):
    """ Returns center scores, removing the flanking ones

    Parameters
    ----------
    all_scores: array_like
        An array containting all scores
    flank_length: int
        Number of flanking sites on each side
    """
    assert flank_length>0
    return all_scores[flank_length:-flank_length]

def get_flanking_scores(all_scores, flank_length):
    """ Returns concatenated flanking scores, removing the center ones

    Parameters
    ----------
    all_scores: array_like
        An array containting all scores
    flank_length: int
        Number of flanking sites on each side
    """
    assert flank_length>0
    return np.concatenate((all_scores[:flank_length], all_scores[-flank_length:]))


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

def create_logo_plot(matplot_dict, meme_dir, logo_path, motif_length):
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

    logo_plot.imshow(logo, extent=[40+15+z*(MAGIC_NUM+1.9),logo.shape[1]+15+XSCALE_FACTOR*(MAGIC_NUM+1.9),0,logo.shape[0]])
    logo_plot.set_axis_off()
    f.add_subplot(logo_plot)
    return logo_plot

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

def create_enrichment_plot(matplot_dict, motif_number, centrimo_txt, centrimo_stats):
    f = matplot_dict['figure']
    gs_h = matplot_dict['gridspec_header']
    gs_b = matplot_dict['gridspec_body']

    all_stats = read_centrimo_stats(centrimo_stats)
    motif_stats = all_stats['MEME']['MOTIF_{}'.format(motif_number)]
    X_values = np.array(motif_stats['pos'])
    print X_values
    Y_values = np.array(motif_stats['count'])
    normalized_Y = Y_values#/np.sum(Y_values)

    centrimo_dict = read_centrimo_txt(centrimo_txt)
    enrichment_pval = float(centrimo_dict['adj_p-value'])
    enrichment = float(centrimo_dict['sites_in_bin'])/float(centrimo_dict['total_sites'])

    enrichment_plot = plt.Subplot(f, gs_h, autoscale_on=True)
    enrichment_plot.set_frame_on(False)
    enrichment_plot.set_xticks([])
    enrichment_plot.set_yticks([])

    enrichment_pval = str('%.1g'%enrichment_pval)
    if 'e' in enrichment_pval:
        enrichment_pval+= '}'
        enrichment_pval= enrichment_pval.replace('e', '*10^{').replace('-0','-')

    textstr = r'\noindent$Enrichment={0:.2f}$\\~\\$(p={1})$'.format(enrichment, enrichment_pval)
    txtx = 0.1*len(textstr)/100.0
    enrichment_plot.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
    f.add_subplot(enrichment_plot)
    enrichment_plot = plt.Subplot(f, gs_b, autoscale_on=True)

    enrichment_plot.plot(X_values, normalized_Y, linewidth=LINEWIDTH)
    enrichment_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')
    enrichment_plot.set_xlabel('$\mathrm{Distance}\ \mathrm{from} \ \mathrm{peak}$', fontsize=FONTSIZE, fontweight='bold')
    enrichment_plot.set_ylabel('$\mathrm{Probability}$', fontsize=FONTSIZE, fontweight='bold')
    enrichment_plot.get_xaxis().tick_bottom()
    enrichment_plot.get_yaxis().tick_left()
    enrichment_plot.get_yaxis().set_tick_params(direction='out')
    enrichment_plot.get_xaxis().set_tick_params(direction='out')
    enrichment_plot.axvline(x=-50, linewidth=3, color='green', linestyle='-.')
    enrichment_plot.axvline(x=50, linewidth=3, color='green', linestyle='-.')

    #enrichment_plot.axvline(x=-100, linewidth=3, color='red', linestyle='-.')
    #enrichment_plot.axvline(x=100, linewidth=3, color='red', linestyle='-.')
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
    for key, cell in table.get_celld().items():
        row, col = key
        if row > 0 and col > 0:
            cell.set_text_props(fontproperties=fontproperties)

    table.set_fontsize(LEGEND_FONTSIZE*8)
    f.add_subplot(ann_plot)

def create_phylop_legend_plot(matplot_dict, motif_freq, sample_phylop_scores, control_phylop_scores, flank_length):
    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']

    phlyop_plots_legend = plt.Subplot(f, gs, autoscale_on=True)
    corr_result = get_pearson_corr(motif_freq, remove_flanking_scores(sample_phylop_scores, flank_length))
    corr_pval = corr_result[1]
    corr_r2 = corr_result[0]
    ttest_result = perform_t_test(remove_flanking_scores(sample_phylop_scores, flank_length), get_flanking_scores(sample_phylop_scores, flank_length))
    p_deltaphylop = ttest_result['one_sided_pval']
    delta_phylop = ttest_result['delta']
    T_deltaphylop = ttest_result['T']
    pearsonr_pval = str('%.1g'%corr_pval)
    if 'e' in pearsonr_pval:
        pearsonr_pval += '}'
        pearsonr_pval = pearsonr_pval.replace('e', '*10^{').replace('-0','-')
    score_pval = str('%.1g'%p_deltaphylop)
    if 'e' in score_pval:
        score_pval += '}'
        score_pval = score_pval.replace('e', '*10^{').replace('-0','-')

    textstr = r'$R_{pearson}=%.2f(p=%s)$' '\n' r'$\Delta_{Phylop}=%.2f(p=%s)$' %(corr_r2, pearsonr_pval, delta_phylop, score_pval)
    txtx = 1-LEGEND_XMULTIPLIER*len(textstr)/100.0
    phlyop_plots_legend.set_frame_on(False)
    phlyop_plots_legend.set_xticks([])
    phlyop_plots_legend.set_yticks([])
    phlyop_plots_legend.text(txtx, TXT_YPOS, textstr, fontsize=LEGEND_FONTSIZE)
    f.add_subplot(phlyop_plots_legend)


def perform_OLS(dependent_variable, independent_variable):
    """Perform Ordinary least square

    Parameters
    ----------

    dependent_variable: array_like
        Array of dependet variable(Y)

    independent_variable: array_like
        Array of indepdent variable(X)

    Returns
    -------

    regression_line: Array like
        Regression line

    """

    regression_fit = sm.OLS(dependent_variable, sm.add_constant(independent_variable)).fit()
    if (len(regression_fit.params)<2):
        ## In cases of multicollinearity simply draw a straight line
        regression_line = independent_variable
    else:
        regression_line = independent_variable*regression_fit.params[1]+regression_fit.params[0]

    return {'regression_fit': regression_fit,
            'regression_line': regression_line}

def create_phylop_scatter(matplot_dict, motif_freq, sample_phylop_scores, control_phylop_scores, flank_length, num_occurrences, y_label):

    f = matplot_dict['figure']
    gs = matplot_dict['gridspec']
    phylop_scatter_plot = plt.Subplot(f, gs, autoscale_on=True)
    control_phylop_scores = remove_flanking_scores(control_phylop_scores, flank_length)
    sample_phylop_scores = remove_flanking_scores(sample_phylop_scores, flank_length)

    fit = np.polyfit(motif_freq, sample_phylop_scores,1)
    fit_fn = np.poly1d(fit)


    control_ols = perform_OLS(control_phylop_scores, motif_freq)
    sample_ols = perform_OLS(sample_phylop_scores, motif_freq)

    sample_regression_line = sample_ols['regression_line']
    control_regression_line = control_ols['regression_line']


    phylop_scatter_plot.scatter(motif_freq, sample_phylop_scores, color='g', s=[POINTSIZE for i in motif_freq])
    phylop_scatter_plot.plot(motif_freq, sample_regression_line, 'g', motif_freq, fit_fn(motif_freq), color='g', linewidth=LINEWIDTH)
    phylop_scatter_plot.scatter(motif_freq, control_phylop_scores, color=GREYNESS, s=[POINTSIZE for i in motif_freq])
    phylop_scatter_plot.plot(motif_freq, control_regression_line, color=GREYNESS, linewidth=LINEWIDTH)

    ticks_and_labels = np.linspace(1.02*min(motif_freq), 1.02*max(motif_freq), num = 5, endpoint=True)
    phylop_scatter_plot.set_xticks(ticks_and_labels)
    ticks_and_labels = ["$%.2f$"%(x/(1.02*num_occurrences)) for x in ticks_and_labels]
    phylop_scatter_plot.set_xticklabels(ticks_and_labels)

    max_yticks = 4
    yloc = plt.MaxNLocator(max_yticks)
    phylop_scatter_plot.yaxis.set_major_locator(yloc)
    phylop_scatter_plot.set_xlabel('$\mathrm{Base}\ \mathrm{Frequency}$', fontsize=FONTSIZE, fontweight='bold')
    phylop_scatter_plot.get_xaxis().tick_bottom()
    phylop_scatter_plot.get_yaxis().tick_left()
    phylop_scatter_plot.set_ylabel('$\mathrm{%s}\ \mathrm{Score}$'%(y_label), fontsize=FONTSIZE, fontweight='bold')
    phylop_scatter_plot.tick_params(axis='y', which='major', pad=TICKPAD)
    phylop_scatter_plot.tick_params(axis='x', which='major', pad=TICKPAD)
    phylop_scatter_plot.get_yaxis().set_tick_params(direction='out')
    phylop_scatter_plot.get_xaxis().set_tick_params(direction='out')
    phylop_scatter_plot.tick_params('both', length=TICKLENGTH, width=2, which='major')

    f.add_subplot(phylop_scatter_plot)


def create_plot(meme_file,
                peak_file,
                fimo_file,
                sample_phylop_file,
                control_phylop_file,
                centrimo_dir,
                sample_gerp_file=None,
                control_gerp_file=None,
                annotate=None,
                motif_number=1,
                flank_length=5):
    meme_record = read_memefile(meme_file)
    record = meme_record['motif_records'][motif_number-1]
    num_occurrences = getattr(record, 'num_occurrences', 'Unknown')
    meme_dir = os.path.abspath(os.path.dirname(meme_file))
    fimo_dir = os.path.abspath(os.path.dirname(fimo_file))

    use_gerp = False
    print sample_gerp_file
    print control_gerp_file

    if annotate == "" or annotate == ' ':
        annotate = None

    if sample_gerp_file:
        use_gerp = True

    max_occur = get_max_occuring_bases(record, max_count=1, count_type=COUNT_TYPE)
    motif_freq = []
    for position in max_occur:
        motif_freq.append(position[0][1])

    motif_freq = np.asarray(motif_freq)
    sample_phylop_scores = np.loadtxt(sample_phylop_file)
    control_phylop_scores = np.loadtxt(control_phylop_file)
    sample_gerp_scores = np.loadtxt(sample_gerp_file)
    control_gerp_scores = np.loadtxt(control_gerp_file)

    motif = record
    motif_length = motif.length
    meme_dir = os.path.abspath(os.path.dirname(meme_file))
    X = [40+15] ## this is by trial and error, the position for the first base logo
    ## Generate all other X coordinates
    for j in range(1,len(motif)+2*flank_length):
        X.append( X[j-1]+MAGIC_NUM+1.9 )

    centrimo_dir = os.path.abspath(centrimo_dir)
    centrimo_txt = os.path.join(centrimo_dir, 'centrimo.txt')
    centrimo_stats = os.path.join(centrimo_dir, 'site_counts.txt')

    ##FIXME This is a big dirty hacl to get thegenerate plots for the Reverse complement logo too
    logo_name =['logo{}.png'.format(motif_number), 'logo_rc{}.png'.format(motif_number)]
    for ln in logo_name:
        setup_matplotlib()
        if 'rc'in ln:
            sample_phylop_scores = sample_phylop_scores[::-1]
        matplot_dict =init_figure(meme_dir=meme_dir, X_values=X, motif=motif_number,
                                    use_gerp=use_gerp, annotate=annotate)
        f = matplot_dict['figure']
        gs = matplot_dict['gs']
        figsize = matplot_dict['figsize']
        right_margin = matplot_dict['right_margin']
        total_px= matplot_dict['total_px']

        logo_plot = create_logo_plot({'figure':f, 'gridspec': gs[0]}, meme_dir, ln, motif_length)


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
        create_stemplot({'figure': f, 'gridspec': gs[1], 'shareX': logo_plot}, X, sample_phylop_scores, motif_length, flank_length=flank_length)


        create_phylop_legend_plot({'figure':f, 'gridspec':gs1[0,0]},  motif_freq, sample_phylop_scores, control_phylop_scores, flank_length)
        create_phylop_scatter({'figure':f, 'gridspec':gs1[1,0]}, motif_freq, sample_phylop_scores, control_phylop_scores, flank_length, num_occurrences, y_label='Phylop')

        if use_gerp:
            create_phylop_legend_plot({'figure':f, 'gridspec':gerp_header_subplot_gs},  motif_freq, sample_gerp_scores, control_gerp_scores, flank_length)
            create_phylop_scatter({'figure':f, 'gridspec':gerp_subplot_gs}, motif_freq, sample_gerp_scores, control_gerp_scores, flank_length, num_occurrences, y_label='Gerp')


        create_enrichment_plot({'figure':f,
                                'gridspec_header':histogram_header_subplot_gs,
                                'gridspec_body': histogram_subplot_gs},
                               motif_number,
                               centrimo_txt,
                               centrimo_stats)

        if 'rc' not in ln:
            out_file = os.path.join(fimo_dir,'motif{}Combined_plots.png'.format(motif_number))
            out_file = 'motif{}Combined_plots.png'.format(motif_number)
        else:
            out_file = os.path.join(fimo_dir,'motif{}Combined_plots_rc.png'.format(motif_number))
            out_file = 'motif{}Combined_plots_rc.png'.format(motif_number)

        if annotate:
            create_annnotation_plot({'figure':f, 'gridspec_header':ann_header_subplot_gs,
                                'gridspec_body': ann_subplot_gs}, annotate)

        f.savefig(out_file, figsize=figsize, dpi=DPI)

if __name__ == '__main__':
    meme_file, peak_file, fimo_file, sample_phylop_file, control_phylop_file, centrimo_dir, sample_gerp_file,  control_gerp_file= sys.argv[1:]
    create_plot(meme_file,
                peak_file,
                fimo_file,
                sample_phylop_file,
                control_phylop_file,
                centrimo_dir,
                sample_gerp_file,
                control_gerp_file,
                annotate=None,
                motif_number=1,
                flank_length=5)
