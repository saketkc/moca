"""Functions to perform stats on data"""
from scipy import stats
from scipy.stats.stats import pearsonr
import numpy as np

def format_pvalue(pval):
    """Latex compatible representations

    """
    pval = str('%.1g' % pval)
    if 'e' in pval:
        pval+= '}'
        pval = pval.replace('e', '*10^{').replace('-0','-')
    return pval

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

def get_pearson_corr(X_values, Y_values):
    """Calculate pearson correlation between inputs

    Arguments
    ---------
    X_values, Y_values: array_like
        Array of values to calculate correaltion on

    Returns
    -------
    pr_p: float
        Pearson correlation
    """
    assert len(X_values) == len(Y_values)
    pr_p = pearsonr(X_values, Y_values)
    return pr_p

def perform_t_test(a, b):
    """Perform independent sample t-test

    Two sided test for the null hypotheisis that two independent
    samples have identical meansi

    $H_0$: Mean score for sample `a` and sample `b` are same
    $H_1$L Mean score for sample `a` > Mean score for sample `b`

    Arguments
    ---------
    a,b: array_like

    Returns
    -------
    statistic: float/array
        Calculate t-statistic

    p-value: float/array
        One sided p-value

    """
    T, two_sided_pval = stats.ttest_ind(a,b)
    delta = np.mean(a)-np.mean(b)
    if T>=0:
        one_sided_pval = two_sided_pval*0.5
    else:
        one_sided_pval = 1-(two_sided_pval*0.5)

    return {'two_sided_pval': two_sided_pval, 'T': T,
            'delta': delta, 'one_sided_pval': one_sided_pval}

def get_center_enrichment(centrimo_txt):
    """Get center encihment and associated p values
    Hypothesis: There is uniform encihment

    """
    centrimo_dict = read_centrimo_txt(centrimo_txt)
    enrichment_pval = float(centrimo_dict['adj_p-value'])
    enrichment = float(centrimo_dict['sites_in_bin'])/float(centrimo_dict['total_sites'])
    return {'enrichment': enrichment, 'enrichment_pval': enrichment_pval}


def perform_regression(X_values, Y_values):
    """Perform linear regression
    Parameters
    ----------
    X_values: array_like

    Y_values: array

    Returns
    -------
    r2: Coefficient of determination

    r2_pval: p-value associated with r2

    """
    pass

def get_motif_evalue(motif):
    """Get motif E-value"""
    return motif.evalue

