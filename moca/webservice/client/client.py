from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from flask import Flask, render_template, request, redirect, url_for
from moca.helpers import ConfigurationParser
from celery import Celery
import sys
import os
from pymongo import MongoClient
from flask import send_from_directory
from celery.result import AsyncResult
from Bio import motifs
import uuid4
from moca.helpers import safe_mkdir
import time
import celery

#TODO organize this

jaspar_motifs = motifs.parse(open('../data/pfm_vertebrates.txt'), 'jaspar')
server_config = None
__form_keys__ = ['assembly']
app = Flask(__name__)

def should_render_json(request):
    req_format = request.args.get('format', '')
    if req_format.lower() == 'json':
        return True
    return False

def api_response(data=None, success=True, message=None):
    return {
                "success": success,
                "data": data,
                "message": message

    }

def insert_task_mapping(task_id, folder_id, metadata):
    collection = db_connector()
    collection.insert_one({'task_id': task_id, 'folder_id': folder_id, 'metadata': metadata})

def db_connector(configuration_file='application.cfg', collection_config_key='mongo_encode_stats_collection'):
    """Connects to MongoDB instance
    Parameters
    ----------
    configuration_file: string
        Path to configuration file

    collection_config_key: string
        Key used in config file for collection name under [mongo] section

    Returns
    -------
    collection: mongodb collection
    """
    config = ConfigurationParser(configuration_file)
    mongo_config = config.get_section('mongo')
    mongo_client = MongoClient('mongodb://{}:{}@{}/{}'.format(mongo_config['mongo_username'],
                                                              mongo_config['mongo_password'],
                                                              mongo_config['mongo_host'],
                                                              mongo_config['mongo_dbname']),
                                    port=int(mongo_config['mongo_port']))
    db = mongo_client.moca_encode_tf
    collection = db[mongo_config['mongo_encode_stats_collection']]
    return collection

def init_params(configuration_file):
    config = ConfigurationParser(configuration_file)
    server_config = config.get_section('server')
    server_config = dict( (key.upper(), value) for key, value in list(server_config.items()))
    for key, value in list(server_config.items()):
        try:
            # For integer
            app.config[key] = int(value)
        except ValueError:
            app.config[key] = value

@app.route('/jaspar')
def jaspar_search(self, tf_name):
    for m in jaspar_motifs:
        if m.name.lower() == tf_name.lower():
            fn = os.path.join(STATIC_PATH, 'logos', m.name+'.png')
            m.weblogo(fn,  show_errorbars=False, logo_title=m.name,  show_fineprint=False , symbols0='A', symbols1='T', symbols2='C', symbols3='G',
                    color0='red', color1='green', color2='blue', color3='orange')
            return api_response(data={'path': m.name+'.png'}, message='png path')
    return api_response()

@app.route('/encoderesults')
def encode_results():
    response = get_all_encode_results()
    data = response.json()
    return render_template('encodeanalysis.html', results=data)

def get_all_encode_results():
    collection = db_connector(collection='encode_data')
    columns_to_return = ['encode_id',
                        'motif_number',
                        'center_enrichment',
                        'delta_phylop_over_control',
                        'delta_phylop_pval',
                        'delta_gerp_over_control',
                        'delta_gerp_pval',
                        'r_phylop_sample',
                        'r_phylop_sample_pval',
                        'r_gerp_sample',
                        'r_gerp_sample_pval',
                        'motif_evalue']
    data = list(collection.find({}, dict((k,1) for k in columns_to_return)))
    for d in data:
        d.update((k, "NaN") for k, v in list(d.items()) if str(v).lower()=='nan')
        d.update((k, float('%.2e' % v)) for k,v in list(d.items()) if 'pval' in k and type(v)==float)
        d.update((k, round(v,3)) for k,v in list(d.items()) if 'pval' not in k and type(v)==float)
        del d['_id']

    return api_response(data=data,
                        message="All records")

def get_encode_metadata(encode_id, motif_number):
    collection = db_connector(collection='encode_data')
    data = list(collection.find({'files.dataset': encode_id}) )
    ## TODO fix this
    assert len(data) == 1
    return api_response(data=data,
                        message="All records")

def get_unique_folder_id():
    return str(uuid4())

def is_request_file_based(request):
    if request.files['file']:
        return True
    return False

def create_job(request):
    folder_id = get_unique_folder_id()
    assembly = request.form['assembly']
    job_folder = os.path.join(server_config['jobs_folder'], folder_id)
    upload_based = is_request_file_based(request)
    safe_mkdir(job_folder)
    if upload_based:
        file_o = request.files['file']
        filename = str(file_o.filename)
        user_filepath = os.path.join(job_folder, filename)
        file_o.save(user_filepath)
    else:
        bed_text = request.form['bedtext']
        filename = 'peaks.bed'
        user_filepath = os.path.join(job_folder, filename)
        with open(user_filepath, 'w') as f:
            f.write(bed_text)
    job_id = submit_job(user_filepath, assembly)
    return job_id

@celery.task(bind=True)
def submit_job(bedfile_path, assembly):
    time.sleep('100')
    return

@app.route('/render/<filename>')
def render_file(filename):
    ##TODO
        return send_from_directory(app.config['UPLOAD_FOLDER'],
                                                                  filename)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'GET':
        return render_template('index.html')
    else:
        job_id = create_job(request)
        return redirect(url_for('results', job_id=job_id))

def get_job_status(job_id):
    ##TODO
    return AsyncResult(job_id).state


def get_job_results(job_id):
    pass

@app.route('/results/<job_id>')
def results(job_id):
    render_json = should_render_json(request)
    status = get_job_status(job_id)
    if status == 'running':
        return api_response(message='running')
    elif status == 'waiting':
        return api_response(message='waiting')
    elif status == 'error':
        return api_response(message='error')
    elif status == 'success' and render_json:
        ## TODO fetch data
        data = get_job_results(job_id)
        return api_response(data=data,message='Done')
    else:
        return render_template('results.html', data=data, should_not_poll=True)

@app.route('/plot/<string:encode_id>/<int:motif_number>')
def get_plot(encode_id, motif_number):
    collection = db_connector(collection='encode_data')
    data = list(collection.find({'encode_id': encode_id,
                                    'motif_number': motif_number}, {'moca_plot': 1, '_id':0}))
    return api_response(data=data,
                        message='Plots')

if __name__ == '__main__':
    init_params(sys.argv[1])
    celery = Celery('MOCA_TASKS',broker=app.config['CELERY_BROKER_URL'])

    #Loads settings for Backend to store results of jobs
    celery.config_from_object('celeryconfig')

    app.run(host='moca.usc.edu', port=8888, debug=True)
