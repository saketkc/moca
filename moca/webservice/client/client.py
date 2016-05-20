from flask import Flask, render_template, request, redirect, url_for, jsonify
from moca.helpers import ConfigurationParser
from celery import Celery
import sys
import requests

app = Flask(__name__)

def init_params(configuration_file):
    config = ConfigurationParser(configuration_file)
    server_config = config.get_section('server')
    server_config = dict( (key.upper(), value) for key, value in server_config.iteritems())
    print server_config
    for key, value in server_config.iteritems():
        try:
            # For integer
            app.config[key] = int(value)
        except ValueError:
            app.config[key] = value

@app.route('/')#encoderesults')
def encode_results():
    response = requests.get('http://moca.usc.edu:8889')
    data = response.json()
    return render_template('encodeanalysis.html', results=data)

if __name__ == '__main__':
    init_params(sys.argv[1])
    celery = Celery('client',
                    broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    app.run(host='moca.usc.edu', port=8888, debug=True)
