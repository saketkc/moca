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
from uuid import uuid4
from moca.helpers import safe_makedir
import time
import celery
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

#TODO organize this

jaspar_motifs = motifs.parse(open('./data/pfm_vertebrates.txt'), 'jaspar')
server_config = None
__form_keys__ = ['assembly']
app = Flask(__name__)

def make_celery(app):
    global server_config
    celery = Celery(app.import_name, broker=server_config['BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task
    class ContextTask(TaskBase):
        abstract = True
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    celery.Task = ContextTask
    return celery


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
    collection = db_connector(collection_name='job_queue')
    collection.insert_one({'task_id': task_id, 'folder_id': folder_id, 'metadata': metadata})

def db_connector(configuration_file='application.cfg', collection_name='encode_tf_stats'):
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
    collection = db[collection_name]
    return collection

def init_params(configuration_file):
    global server_config
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
    collection = db_connector(collection='encode_data', collection_name='encode_tf_stats')
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
    collection = db_connector(collection_name='tf_metadata')
    data = list(collection.find({'files.dataset': encode_id}) )
    ## TODO fix this
    assert len(data) == 1
    return api_response(data=data,
                        message="All records")

def get_unique_folder_id():
    return str(uuid4())

def is_request_file_based(request):
    print(request.files)
    if request.files['bedfile']:
        return True
    return False

def create_job(request):
    global server_config
    print(request.form)
    folder_id = get_unique_folder_id()
    print (folder_id)
    assembly = request.form['assembly']
    print (assembly)
    print(server_config)
    job_folder = os.path.join(server_config['JOBS_FOLDER'], folder_id)
    print (job_folder)
    upload_based = is_request_file_based(request)
    print(upload_based)
    safe_makedir(job_folder)
    print('done')
    if upload_based:
        file_o = request.files['bedfile']
        filename = str(file_o.filename)
        user_filepath = os.path.join(job_folder, filename)
        file_o.save(user_filepath)
    else:
        bed_text = request.form['bedtext']
        filename = 'peaks.bed'
        user_filepath = os.path.join(job_folder, filename)
        with open(user_filepath, 'w') as f:
            f.write(bed_text)
    job_id = submit_job.apply_async([user_filepath, assembly])
    return job_id

@celery.task(name='client.submit_job')
def submit_job(self, data):
    time.sleep(1000)
    return True

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

@app.route('/results/<job_id>')
def results(job_id):
    render_json = should_render_json(request)
    status = get_job_status(job_id)
    print('STATUS: {}'.format(status))
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
        data = []
        return render_template('results.html', data=data, should_not_poll=True)

@app.route('/plot/<string:encode_id>/<int:motif_number>')
def get_plot(encode_id, motif_number):
    collection = db_connector(collection_name='encode_data')
    data = list(collection.find({'encode_id': encode_id,
                                 'motif_number': motif_number},
                                {'moca_plot': 1, '_id':0}))
    return api_response(data=data,
                        message='Plots')

if __name__ == '__main__':
    init_params(sys.argv[1])
    #celery = Celery('client',
    #                broker='mongodb://localhost:27017/celery_broker',
    #                backend='mongodb')

    CELERY_IMPORTS=('client',)

    CELERY_RESULT_BACKEND = "mongodb"
    CELERY_MONGODB_BACKEND_SETTINGS = {
        "host": "127.0.0.1",
        "port": 27017,
        "database": "celery_jobs",
        "taskmeta_collection": "stock_taskmeta_collection",
    }
   # celery.config_from_object('celeryconfig')
    config = {}
    config['BROKER_URL'] = server_config['BROKER_URL']
    config['CELERY_RESULT_BACKEND'] = 'mongodb'
    config['CELERY_MONGODB_BACKEND_SETTINGS'] = CELERY_MONGODB_BACKEND_SETTINGS
    config['CELERY_IMPORTS'] = ('client', )

    #celery.conf.update(config)
    app.config.update(config)
    celery = make_celery(app)



    app.run(host='moca.usc.edu', port=8888, debug=True)
