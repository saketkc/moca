from __future__ import division
import os
from pymongo import MongoClient
import matplotlib.pyplot as plt
import searborn

__root_dir__ = '/media/data1/encode_analysis'
flank_length = 5
COUNT_TYPE = 'counts'

client = MongoClient()
db = client.moca_encode_tf
collection = db.encode_tf_stats
cursor = collection.find()

for documentcursor = collection.find()c

plt.hist(collection.find())
