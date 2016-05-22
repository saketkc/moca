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

class Webservice(Resource):
    """Base class to implement API"""
    def __init__(self, configuration_file):
        """Init

        Parameters
        ----------
        configuration_file: str
            Path to configuration file
        """
        self.config = ConfigurationParser(configuration_file)
        self.mongo_config = self.config.get_section('mongo')
        self.mongo_client = MongoClient('mongodb://{}:{}@{}/{}'.format(self.mongo_config['mongo_username'],
                                                                    self.mongo_config['mongo_password'],
                                                                    self.mongo_config['mongo_host'],
                                                                    self.mongo_config['mongo_dbname']),
                                        port=int(self.mongo_config['mongo_port']))
        self.db = self.mongo_client.moca_encode_tf
        self.collection = self.db[self.mongo_config['mongo_collection']]

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
                                  'motif_evalue',
                                  'motif_logo',
                                  'motif_logorc',
                                  'moca_plot',
                                  'moca_plotrc']

    def get(self):
        """Get all records from collection"""
        data = list(self.collection.find())
        for d in data:
            d.update((k, "NaN") for k, v in d.iteritems() if str(v).lower()=='nan')
            temp_d = {}
            for k in self.columns_to_return:
                if k not in d:
                    temp_d[k] = 'NAN'
                else:
                    temp_d[k]= d[k]
            keys = d.keys()
            for key in keys:
                if key not in self.columns_to_return:
                    del d[key]
            d.update(temp_d)
            d.update((k, float('%.2e' % v)) for k,v in d.iteritems() if 'pval' in k and type(v)==float)
            d.update((k, round(v,3)) for k,v in d.iteritems() if 'pval' not in k and type(v)==float)

            #for key in d.keys():
            #    if key in ['gerp_mean_sample', 'gerp_mean_control', 'phylop_mean_sample', 'phylop_mean_control']:
            #        del d[key]

        print data[0]
        aoColumns = []
        for key in self.columns_to_return:
            aoColumns.append({"mDataProp": key, "sTitle": key})

        return api_response(data={'rows': data, 'columns': aoColumns},
                            success="success",
                            message="All records")

    def get_filtered_records(self, filter_dict=None):
        """Get records post filtering

        Parameters
        ----------

        filter_dict: dict
            A dictionary with key=filter key,

        """
        pass



if __name__ == '__main__':
    if len(sys.argv)!=2:
        print  'Run: python api.py <configuration.cfg>'
        sys.exit(1)
    print sys.argv[1]
    app = Flask(__name__)
    api = Api(app)
    api.add_resource(Webservice, '/', resource_class_args = [sys.argv[1]])
    app.run(host='moca.usc.edu',
            debug=True,
            port=8889)
