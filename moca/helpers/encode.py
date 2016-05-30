#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import os
import requests
import shutil
from pymongo import MongoClient
from moca.helpers import safe_makedir
import io
import gzip
import json

__base_url__ = 'https://www.encodeproject.org/'
__root_dir__ = '/media/data1/ENCODE_V3/'
ALLOWED_OUTPUT_TYPES = ['optimal idr thresholded peaks', 'peaks']
ALLOWED_FILETYPES = ['bed narrowPeak', 'bed broadPeak']

def download_peakfile(source_url, filename, destination_dir):
    """Download peakfile from encode"""
    response = requests.get(source_url, stream=True)
    with open(os.path.join(destination_dir, filename), 'wb') as f:
        shutil.copyfileobj(response.raw, f)

    with gzip.open(os.path.join(destination_dir, filename), 'rb') as in_file:
        with open(os.path.join(destination_dir, filename.replace('.gz','')), 'wb') as out_file:
            out_file.write( in_file.read()  )
    del response

def download_idr_tfs(root_dir, metadata):
    """Download all tfs with idr called peaks"""
    idr_records = fetch_idr_record(metadata)
    ## Theere is only one IDR per sample
    if len(idr_records)!=1:
        print(idr_records[0]['dataset'])
    assert len(idr_records) <= 1
    for idr_record in idr_records:
        dataset = idr_record['dataset']
        peakfilename = idr_record['peakfilename'] + '.bed.gz'
        dataset_dir = os.path.join(root_dir, dataset)
        safe_makedir(dataset_dir)
        source_url = __base_url__ + idr_record['href']
        print(source_url)
        download_peakfile(source_url, peakfilename, dataset_dir)
        save_metadata_json(idr_record, dataset_dir)
        return {'assembly': idr_record['assembly'],'bedfile': os.path.join(dataset_dir, peakfilename.replace('.gz',''))}

def save_metadata(metadata):
    """Save metadata to mongodb"""
    client = MongoClient()
    db = client.moca_encode_tf
    result = db.tf_metadata.insert_one(metadata)

def save_metadata_json(metadata, directory):
    """Save metadata locally"""
    with open(os.path.join(directory,'metadata.json'), 'w') as outfile:
        json.dump(metadata, outfile)

def get_experiment(experiment_id):
    """Get and save metadata for an experiment"""
    req = requests.get("{}experiments/{}/?format=json".format(__base_url__, experiment_id))
    metadata = req.json()
    if metadata['status'] != 'error':
        download_idr_tfs(__root_dir__, metadata)
    else:
        print('Error fetching metadata for {}'.format(experiment_id))

def get_idr_controlled_peaks():
    """Return records that have atleast one idr called peak"""
    client = MongoClient()
    db = client.moca_encode_tf
    results = db.tf_metadata_0529016.find({'files.output_type': 'optimal idr thresholded peaks', 'assembly': 'hg19'}, no_cursor_timeout=True)
    data = results[:]
    client.close()
    return data

def fetch_idr_record(metadata):
    files = metadata['files']
    biosample_term_name = metadata['biosample_term_name']
    assay_term_name = metadata['assay_term_name']
    description = metadata['description']
    gene_name = metadata['target']['label']
    parent_metadata = {'biosample_term_name': biosample_term_name,
                       'assay_term_name': assay_term_name,
                       'description': description,
                       'gene_name': gene_name}
    idr_records = []
    for f in files:
        file_status = f['status']
        file_type = f['file_type']
        output_type = f['output_type']
        if output_type == 'optimal idr thresholded peaks' and file_type in ALLOWED_FILETYPES and file_status == 'released':
            dataset = f['dataset']
            dataset = dataset.replace('experiments','').replace('/','')
            href = f['href']
            title = f['title']
            assembly = f['assembly']
            idr_records.append({'href': href, 'metadata':f, 'parent_metadata': parent_metadata, 'dataset': dataset, 'peakfilename': title, 'assembly': assembly})
    return idr_records

def get_encode_peakfiles(encode_id):
    req = requests.get("{}experiments/{}/?format=json".format(__base_url__, encode_id))
    metadata = req.json()
    status = metadata['status']
    if status == 'error':
        return metadata
    files = metadata['files']
    biosample_term_name = metadata['biosample_term_name']
    assay_term_name = metadata['assay_term_name']
    description = metadata['description']
    gene_name = metadata['target']['label']
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
    req = requests.get("{}experiments/{}/?format=json".format(__base_url__, dataset))
    metadata = req.json()
    status = metadata['status']
    files = metadata['files']
    biosample_term_name = metadata['biosample_term_name']
    assay_term_name = metadata['assay_term_name']
    description = metadata['description']
    gene_name = metadata['target']['label']
    if status == 'error':
        return metadata
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
        print(experiment)
        get_experiment(experiment)

if __name__ == '__main__':
    list_of_metadata = get_idr_controlled_peaks()
    for data in list(list_of_metadata):
        download_idr_tfs(__root_dir__, data)

