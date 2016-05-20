import sys
from pymongo import MongoClient
from flask import Flask
from flask_restful import Resource, Api
from flask import render_template
from flask_restful import reqparse

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

    def get(self):
        """Get all records from collection"""
        data = list(self.collection.find())
        map(lambda d: d.pop('_id'), data)
        for d in data:
            d.update((k, "NaN") for k, v in d.iteritems() if str(v).lower()=='nan')
            for key in d.keys():
                if key in ['gerp_mean_sample', 'gerp_mean_control', 'phylop_mean_sample', 'phylop_mean_control']:
                    del d[key]

        print data[0]
        all_column_names = []
        for d in data:
            for key in d.keys():
                if key not in all_column_names:
                    all_column_names.append(key)

        aoColumns = []
        for key in all_column_names:
            aoColumns.append({"mDataProp": key, "sTitle": key})
        print all_column_names
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
