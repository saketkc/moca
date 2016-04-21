from __future__ import division
import os
from pymongo import MongoClient
import numpy as np
from moca.helpers import read_memefile, get_total_sequences
from moca.helpers import create_binary_pickle, unpickle_numpy_array
from moca.helper.seqstats import get_center_enrichment, get_motif_evalue, perform_t_test, get_pearson_corr, get_flanking_sequences, remove_flanking_sequences

__root_dir__ = '/media/data1/encode_analysis'

for d in os.listdir(__root_dir__):

    client = MongoClient()
    db = client.moca_encode_tf
    results = db.tf_metadata.find({'@id': '/experiments/{}/'.format(d)})
    meme_file = os.path.join(__root_dir__, d, 'moca_output', 'meme_out', 'meme.txt')
    if not os.path.isfile(meme_file):
        print 'Skipping {}'.format(d)
        continue
    meme_info = read_memefile(meme_file)
    sample_phylop_scores = np.loadtxt(sample_phylop_file)
    control_phylop_scores = np.loadtxt(control_phylop_file)
    sample_gerp_scores = np.loadtxt(sample_gerp_file)
    control_gerp_scores = np.loadtxt(control_gerp_file)
    total_sequences = get_total_sequences(meme_file)
    for i in range(0, meme_info['total_motifs']):
        fimo_main = os.path.join(os.path.dirname(meme_file), 'fimo_out_{}'.format(i+1))
        fimo_random = os.path.join(os.path.dirname(meme_file), 'fimo_random_{}'.format(i+1))

        gerp_mean_sample = np.loadtxt(os.path.join(fimo_main, 'gerp.mean.txt')).tolist()
        phylop_mean_sample = np.loadtxt(os.path.join(fimo_main, 'phylop.mean.txt')).tolist()
        gerp_mean_random= np.loadtxt(os.path.join(fimo_random, 'gerp.mean.txt')).tolist()
        phylop_mean_random = np.loadtxt(os.path.join(fimo_random, 'phylop.mean.txt')).tolist()

        gerp_raw= np.loadtxt(os.path.join(fimo_main, 'gerp.raw.txt'))
        phylop_raw = np.loadtxt(os.path.join(fimo_main, 'phylop.raw.txt'))

        motif_enrichment = meme_info['motif_occurrences']['motif{}'.format(i+1)]/total_sequences

        delta_phylop_ttest =
    ttest_result = perform_t_test(remove_flanking_scores(sample_phylop_scores, flank_length), get_flanking_scores(sample_phylop_scores, flank_length))
    p_deltaphylop = ttest_result['one_sided_pval']
    delta_phylop = ttest_result['delta']
        db.encode_tf_stats.insert_one({ 'encode_id': d,
                                        'motif_number': i+1,
                                        'gerp_mean': gerp_mean,
                                        'phylop_mean': phylop_mean,
                                        'gerp_raw': create_binary_pickle(gerp_raw),
                                        'phylop_raw': create_binary_pickle(phylop_raw),
                                        'motif_enrichment': motif_enrichment,
                                        'delta_phylop_over_control': delta_phylop,
                                        'delta_gerp_over_control': delta_gerp,
                                        'r2_phylop_main': r2_phylop_main,
                                        'r2_phylop_main_pval': r2_phylop_main_pval,
                                        'r2_phylop_control': r2_phylop_control,
                                        'r2_phylop_control_pval': r2_phylop_control_pval,
                                        'r2_gerp_main': r2_gerp_main,
                                        'r2_gerp_control': r2_gerp_control,
                                        'center_enrichment': center_enrichment,
                                        'center_enrichment_pval': center_enrichment_pval,
                                        'motif_evalue': motif_evalue
                                    })

