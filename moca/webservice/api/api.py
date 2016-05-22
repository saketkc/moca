import sys
from pymongo import MongoClient
from flask import Flask
from flask_restful import Resource, Api

from moca.helpers import ConfigurationParser

def api_response(data, success, message):
    return {
                "success": success,
                "data": data,
                "message": message

    }

#TODO Fix this and the method below

def db_connector(configuration_file):
    config = ConfigurationParser(configuration_file)
    mongo_config = config.get_section('mongo')
    mongo_client = MongoClient('mongodb://{}:{}@{}/{}'.format(mongo_config['mongo_username'],
                                                              mongo_config['mongo_password'],
                                                              mongo_config['mongo_host'],
                                                              mongo_config['mongo_dbname']),
                                    port=int(mongo_config['mongo_port']))
    db = mongo_client.moca_encode_tf
    collection = db[mongo_config['mongo_collection']]
    return collection

class Webservice(Resource):
    """Base class to implement API"""
    def __init__(self, configuration_file):
        """Init

        Parameters
        ----------
        configuration_file: str
            Path to configuration file
        """
        self.collection = db_connector(configuration_file)

        self.columns_to_return = ['encode_id',
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

    def get(self):
        """Get all records from collection"""
        data = list(self.collection.find({}, dict((k,1) for k in self.columns_to_return)))
        for d in data:
            d.update((k, "NaN") for k, v in d.iteritems() if str(v).lower()=='nan')
            d.update((k, float('%.2e' % v)) for k,v in d.iteritems() if 'pval' in k and type(v)==float)
            d.update((k, round(v,3)) for k,v in d.iteritems() if 'pval' not in k and type(v)==float)
            del d['_id']
        aoColumns = []
        for key in self.columns_to_return:
            aoColumns.append({"mDataProp": key, "sTitle": key})

        return api_response(data={'rows': data, 'columns': aoColumns},
                            success="success",
                            message="All records")


class GetPlot(Resource):
    def __init__(self, configuration_file):
        """Init

        Parameters
        ----------
        configuration_file: str
            Path to configuration file
        """
        self.collection = db_connector(configuration_file)

    def get(self, encode_id, motif_number):
        data = list(self.collection.find({'encode_id': encode_id,
                                          'motif_number': motif_number}, {'moca_plot': 1, '_id':0}))
        return api_response(data={'rows': data, 'columns': {}},
                            success="success",
                            message="All records")


if __name__ == '__main__':
    if len(sys.argv)!=2:
        print  'Run: python api.py <configuration.cfg>'
        sys.exit(1)
    app = Flask(__name__)
    api = Api(app)
    api.add_resource(Webservice, '/', resource_class_args = [sys.argv[1]])
    api.add_resource(GetPlot, '/plot/<string:encode_id>/<int:motif_number>', resource_class_args = [sys.argv[1]])
    app.run(host='moca.usc.edu',
            debug=True,
            port=8889)
