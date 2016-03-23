#!/usr/bin/env python

import os
import requests
import shutil
from pymongo import MongoClient

__base_url__ = 'https://www.encodeproject.org/'

ALLOWED_OUTPUT_TYPES = ['optimal idr thresholded peaks', 'peaks']
ALLOWED_FILETYPES = ['bed narrowPeak', 'bed broadPeak']

def download_peakfile(source_url, filename, destination_dir):
    response = requests.get(source_url, stream=True)
    with open(os.path.join(destination_dir, filename), 'wb') as f:
        shutil.copyfileobj(response.raw, f)
    del response

def save_metadata(metadata):
    client = MongoClient()
    db = client.moca_encode_tf
    result = db.tf_metadata.insert_one(metadata)
    print result.inserted_id

def get_experiment(experiment_id):
    req = requests.get("{}experiments/{}/?format=json".format(__base_url__, experiment_id))
    response_json = req.json()
    if response_json['status'] != 'error':
        save_metadata(response_json)
    else:
        print 'Error fetching metadata for {}'.format(experiment_id)


def get_idr_controlled_peaks():
    client = MongoClient()
    db = client.moca_encode_tf
    results = db.tf_metadata.find({'files.output_type': 'optimal idr thresholded peaks'})
    return results

def get_encode_peakfiles(encode_id):
    req = requests.get("{}experiments/{}/?format=json".format(__base_url__, encode_id))
    response_json = req.json()
    status = response_json['status']
    if status == 'error':
        return response_json
    files = response_json['files']
    biosample_term_name = response_json['biosample_term_name']
    assay_term_name = response_json['assay_term_name']
    description = response_json['description']
    gene_name = response_json['target']['label']
    files_to_download = []
    for f in files:
        file_type = f['file_type']
        file_status = f['status']
        if file_type in ALLOWED_FILETYPES:
            dataset = f['dataset']
            dataset = dataset.replace('experiments','').replace('/','')
            assert (dataset == encode_id)
            assembly = f['assembly']
            href = f['href']
            output_type = f['output_type']
            try:
               tech_repl_number =  f['replicate']['technical_replicate_number']
               bio_repl_number = f['replicate']['biological_replicate_number']
            except KeyError:
                assert (output_type in ALLOWED_OUTPUT_TYPES)
                tech_repl_number =  ''
                bio_repl_number = ''

            files_to_download.append({'href': href, 'tech_repl_number': tech_repl_number,
                                      'bio_repl_number': bio_repl_number,
                                      'output_type': output_type.replace('bed ',''),
                                      'file_type': file_type,
                                      'title': f['title'],
                                      'dataset':dataset,
                                      'assembly': assembly,
                                      'biosample_term_name': biosample_term_name,
                                      'assay_term_name': assay_term_name,
                                      'gene_name': gene_name,
                                      'description': description,
                                      'file_status': file_status})
    return files_to_download

def get_metadata_for_peakfile(dataset, peakfile):
    #TODO this is repeated from last function and probaly is an inefficient way to implement this
    req = requests.get("{}experiments/{}/?format=json".format(__base_url__, dataset))
    response_json = req.json()
    status = response_json['status']
    files = response_json['files']
    biosample_term_name = response_json['biosample_term_name']
    assay_term_name = response_json['assay_term_name']
    description = response_json['description']
    gene_name = response_json['target']['label']
    if status == 'error':
        return response_json
    for f in files:
        title = f['title']
        file_status = f['status']
        if title == peakfile:
            assembly = f['assembly']
            href = __base_url__ + f['href']
            output_type = f['output_type']
            try:
               tech_repl_number =  f['replicate']['technical_replicate_number']
               bio_repl_number = f['replicate']['biological_replicate_number']
            except KeyError:
                assert (output_type in ALLOWED_OUTPUT_TYPES)
                tech_repl_number =  ''
                bio_repl_number = ''
            return {'href': href,
                    'tech_repl_number': tech_repl_number,
                    'bio_repl_number': bio_repl_number,
                    'output_type': output_type,
                    'title': f['title'],
                    'file_type': f['file_type'].replace('bed ',''),
                    'dataset':dataset,
                    'assembly': assembly,
                    'biosample_term_name': biosample_term_name,
                    'assay_term_name': assay_term_name,
                    'gene_name': gene_name,
                    'description': description,
                    'file_status': file_status,
                    }

    return None

def search_encode_tfs():
    url = __base_url__ + '/search?type=Experiment&assay_title=ChIP-seq&limit=all&status=released&target.investigated_as=transcription+factor&format=json'
    req = requests.get(url)
    resp_json = req.json()
    all_samples = resp_json['@graph']
    all_tfs = [sample['target']['label'] for sample in all_samples]
    all_tfs =  set(all_tfs)
    all_experiments = [sample['@id'].strip().replace('/','').replace('experiments', '') for sample in all_samples]
    for experiment in all_experiments:
        print experiment
        get_experiment(experiment)

