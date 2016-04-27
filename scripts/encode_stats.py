from __future__ import division
import os
from pymongo import MongoClient
import numpy as np
from moca.helpers import read_memefile, get_total_sequences
from moca.helpers.db import create_binary_pickle#, unpickle_numpy_array
from moca.helpers import get_max_occuring_bases
#from moca.helpers import read_centrimo_txt
#from moca.helpers import read_centrimo_stats
from moca.helpers.seqstats import get_center_enrichment, get_motif_evalue, perform_t_test, get_pearson_corr
from moca.helpers.seqstats import get_flanking_scores, remove_flanking_scores, perform_OLS

__root_dir__ = '/media/data1/encode_analysis'
flank_length = 5
COUNT_TYPE = 'counts'

client = MongoClient()
db = client.moca_encode_tf
db.encode_tf_stats.remove()

for d in os.listdir(__root_dir__):

    results = db.tf_metadata.find({'@id': '/experiments/{}/'.format(d)})
    meme_file = os.path.join(__root_dir__, d, 'moca_output', 'meme_out', 'meme.txt')
    centrimo_dir = os.path.join(__root_dir__, d, 'moca_output', 'centrimo_out')
    if not os.path.isfile(meme_file):
        print 'Skipping {}'.format(d)
        continue
    meme_info = read_memefile(meme_file)

    total_sequences = get_total_sequences(meme_file)

    for i in range(0, meme_info['total_motifs']):
        record = meme_info['motif_records'][i]
        max_occur = get_max_occuring_bases(record, max_count=1, count_type=COUNT_TYPE)
        motif_freq = []
        for position in max_occur:
            motif_freq.append(position[0][1])

        motif_freq = np.asarray(motif_freq)

        fimo_sample = os.path.join(os.path.dirname(meme_file), 'fimo_out_{}'.format(i+1))
        fimo_random = os.path.join(os.path.dirname(meme_file), 'fimo_random_{}'.format(i+1))

        if not os.path.isfile(os.path.join(fimo_sample, 'gerp.mean.txt')):
            db.encode_tf_stats.insert_one({ 'encode_id': d,
                                            'motif_number': i+1,
                                            'motif_missing_error': True})
            continue

        motif_enrichment = meme_info['motif_occurrences']['motif{}'.format(i+1)]/total_sequences
        centrimo_dir = os.path.abspath(centrimo_dir)
        centrimo_txt = os.path.join(centrimo_dir, 'centrimo.txt')
        enrichment_info = get_center_enrichment(centrimo_txt, i+1)
        center_enrichment = enrichment_info['enrichment']
        center_enrichment_pval = enrichment_info['enrichment_pval']
        motif_evalue = get_motif_evalue(record)

        if os.stat(os.path.join(fimo_sample, 'gerp.mean.txt')).st_size == 0:
            db.encode_tf_stats.insert_one({ 'encode_id': d,
                                            'motif_number': i+1,
                                            'center_enrichment': center_enrichment,
                                            'center_enrichment_pval': center_enrichment_pval,
                                            'motif_evalue': motif_evalue,
                                            'motif_enrichment': motif_enrichment,
                                            'no_fimo_hit_sample': True})
            continue

        gerp_mean_sample = np.loadtxt(os.path.join(fimo_sample, 'gerp.mean.txt')).tolist()
        phylop_mean_sample = np.loadtxt(os.path.join(fimo_sample, 'phylop.mean.txt')).tolist()

        delta_phylop_ttest = perform_t_test(remove_flanking_scores(phylop_mean_sample, flank_length),
                                            get_flanking_scores(phylop_mean_sample, flank_length))
        p_delta_phylop = delta_phylop_ttest['one_sided_pval']
        delta_phylop = delta_phylop_ttest['delta']

        delta_gerp_ttest = perform_t_test(remove_flanking_scores(gerp_mean_sample, flank_length),
                                          get_flanking_scores(gerp_mean_sample, flank_length))
        p_delta_gerp = delta_gerp_ttest['one_sided_pval']
        delta_gerp = delta_gerp_ttest['delta']

        phylop_sample_ols = perform_OLS(remove_flanking_scores(phylop_mean_sample, flank_length), motif_freq)
        gerp_sample_ols = perform_OLS(remove_flanking_scores(gerp_mean_sample, flank_length), motif_freq)
        phylop_sample_fit = phylop_sample_ols['regression_fit']
        gerp_sample_fit = gerp_sample_ols['regression_fit']
        corr_phylop_sample = get_pearson_corr(motif_freq, remove_flanking_scores(phylop_mean_sample, flank_length))
        corr_gerp_sample = get_pearson_corr(motif_freq, remove_flanking_scores(gerp_mean_sample, flank_length))
        r_phylop_sample, r_phylop_sample_pval = corr_phylop_sample
        r_gerp_sample, r_gerp_sample_pval = corr_gerp_sample
        r2_phylop_sample = phylop_sample_fit.rsquared
        r2_gerp_sample = gerp_sample_fit.rsquared


        #gerp_raw= np.loadtxt(os.path.join(fimo_sample, 'gerp.raw.txt'))
        #phylop_raw = np.loadtxt(os.path.join(fimo_sample, 'phylop.raw.txt'))

        if os.stat(os.path.join(fimo_random, 'gerp.mean.txt')).st_size == 0:
            db.encode_tf_stats.insert_one({ 'encode_id': d,
                                            'motif_number': i+1,
                                            'center_enrichment': center_enrichment,
                                            'center_enrichment_pval': center_enrichment_pval,
                                            'motif_evalue': motif_evalue,
                                            'motif_enrichment': motif_enrichment,
                                            'phylop_mean_sample': phylop_mean_sample,
                                            'gerp_mean_sample': gerp_mean_sample,
                                            'r2_phylop_sample': r2_phylop_sample,
                                            'r_phylop_sample': r_phylop_sample,
                                            'r_phylop_sample_pval': r_phylop_sample_pval,
                                            'r2_gerp_sample': r2_gerp_sample,
                                            'r_gerp_sample': r_gerp_sample,
                                            'r_gerp_sample_pval': r_gerp_sample_pval,
                                            'no_fimo_hit_control': True})
            continue


        gerp_mean_control =  np.loadtxt(os.path.join(fimo_random, 'gerp.mean.txt')).tolist()
        phylop_mean_control = np.loadtxt(os.path.join(fimo_random, 'phylop.mean.txt')).tolist()
        phylop_control_ols = perform_OLS(remove_flanking_scores(phylop_mean_control, flank_length), motif_freq)
        gerp_control_ols = perform_OLS(remove_flanking_scores(gerp_mean_control, flank_length), motif_freq)
        phylop_control_fit = phylop_control_ols['regression_fit']
        gerp_control_fit = gerp_control_ols['regression_fit']
        corr_phylop_control = get_pearson_corr(motif_freq, remove_flanking_scores(phylop_mean_control, flank_length))
        corr_gerp_control = get_pearson_corr(motif_freq, remove_flanking_scores(gerp_mean_control, flank_length))
        r_phylop_control, r_phylop_control_pval = corr_phylop_control
        r_gerp_control, r_gerp_control_pval = corr_gerp_control
        r2_phylop_control = phylop_control_fit.rsquared
        r2_gerp_control = gerp_control_fit.rsquared



        db.encode_tf_stats.insert_one({ 'encode_id': d,
                                        'motif_number': i+1,
                                        'center_enrichment': center_enrichment,
                                        'center_enrichment_pval': center_enrichment_pval,
                                        'motif_evalue': motif_evalue,
                                        'motif_enrichment': motif_enrichment,
                                        'phylop_mean_sample': phylop_mean_sample,
                                        'gerp_mean_sample': gerp_mean_sample,
                                        'r2_phylop_sample': r2_phylop_sample,
                                        'r_phylop_sample': r_phylop_sample,
                                        'r_phylop_sample_pval': r_phylop_sample_pval,
                                        'r2_gerp_sample': r2_gerp_sample,
                                        'r_gerp_sample': r_gerp_sample,
                                        'r_gerp_sample_pval': r_gerp_sample_pval,

                                        #'gerp_raw': create_binary_pickle(gerp_raw),
                                        #'phylop_raw': create_binary_pickle(phylop_raw),
                                        'phylop_mean_control': phylop_mean_control,
                                        'gerp_mean_control': gerp_mean_control,
                                        'delta_phylop_over_control': delta_phylop,
                                        'delta_phylop_pval': p_delta_phylop,
                                        'delta_gerp_over_control': delta_gerp,
                                        'delta_gerp_pval': p_delta_gerp,
                                        'r2_phylop_control': r2_phylop_control,
                                        'r_phylop_control': r_phylop_control,
                                        'r_phylop_control_pval': r_phylop_control_pval,
                                        'r2_gerp_control': r2_gerp_control,
                                        'r_gerp_control_pval': r_gerp_control_pval,
                                    })

