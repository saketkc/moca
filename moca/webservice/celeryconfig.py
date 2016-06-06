#from celery.schedules import crontab

CELERY_IMPORTS=('client', )
CELERY_BROKER_URL = 'mongodb://localhost:27017/celery_broker'
BROKER_URL = 'mongodb://localhost:27017/celery_broker'

CELERY_RESULT_BACKEND = "mongodb"
CELERY_MONGODB_BACKEND_SETTINGS = {
    "host": "127.0.0.1",
    "port": 27017,
    "database": "celery_jobs",
    "taskmeta_collection": "stock_taskmeta_collection",
}

#used to schedule tasks periodically and passing optional arguments
#Can be very useful. Celery does not seem to support scheduled task but only periodic
#CELERYBEAT_SCHEDULE = {
#    'every-minute': {
#        'task': 'client.submit_job',
#        'schedule': crontab(minute='*/1'),
#        'args': (1,2),
#    },
#}
