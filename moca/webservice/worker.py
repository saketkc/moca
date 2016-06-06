from __future__ import absolute_import, unicode_literals

from celery import current_app
from celery.bin import worker


if __name__ == '__main__':
    app = current_app._get_current_object()
    options = {
        'broker': 'mongodb://localhost:27017/celery_broker',
        'loglevel': 'INFO',
        'traceback': True,
    }

    worker = worker.worker(app=app)#, broker=options['broker'])


    worker.run(**options)

