import os
from pymongo import MongoClient

__root_dir__ = '/media/data1/encode_analysis'

for d in os.listdir(__root_dir__):

    client = MongoClient()
    db = client.moca_encode_tf
    results = db.tf_metadata.find({'@id': '/experiments/{}/'.format(d)})
    print d
    print results.count()
