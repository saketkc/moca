from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import sys
from eve import Eve
from moca.helpers import ConfigurationParser

def get_startup_params(config_file):
    parsed_config = ConfigurationParser(config_file)
    return parsed_config.get_section('mongo')

def run_server(configuration_file):
    settings = get_startup_params(configuration_file)
    settings = dict((k.upper(), v) for k, v in settings.items())
    print(settings)
    settings['DOMAIN'] = {'encode_tf_stats': {}}# {'datasource': {'source': 'moca_encode_tf'}}}
    settings['ALLOW_UNKNOWN'] = True
    settings['DEBUG'] = True
    app = Eve(settings=settings)
    app.run(host='moca.usc.edu', port=8888)

if __name__ == '__main__':
    run_server(sys.argv[1])



