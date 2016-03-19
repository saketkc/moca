"""Functions to perform stats on data"""
from scipy import stats
from scipy.stats.stats import pearsonr
import numpy as np

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

