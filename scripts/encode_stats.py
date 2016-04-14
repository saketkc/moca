from __future__ import division
import os
from pymongo import MongoClient
import numpy as np
from moca.helpers import read_memefile, get_total_sequences
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
    total_sequences = get_total_sequences(meme_file)
    for i in range(0, meme_info['total_motifs']):
        fimo_main = os.path.join(os.path.dirname(meme_file), 'fimo_out_{}'.format(i+1))
        fimo_random = os.path.join(os.path.dirname(meme_file), 'fimo_random_{}'.format(i+1))
        motif_enrichment = meme_info['motif_occurrences']['motif{}'.format(i+1)]/total_sequences

        gerp_mean = np.loadtxt(os.path.join(fimo_main, 'gerp.mean.txt'))
        phylop_mean = np.loadtxt(os.path.join(fimo_main, 'phylop.mean.txt'))

        if motif_enrichment>0.5 and phylop_mean > 1.5 and
            db.good_tf.insert_one({''})
        else:

