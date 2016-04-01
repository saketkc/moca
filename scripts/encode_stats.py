import os
import pymongo

__root_dir__ = '/media/data1/encode_analysis'

for d in os.listdir(__root_dir__):

    client = MongoClient()
    db = client.moca_encode_tf
    results = db.tf_metadata.find({'files.output_type': 'optimal idr thresholded peaks'}, no_cursor_timeout=True)
    data = results[:]
    client.close()
