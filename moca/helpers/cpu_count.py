from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import multiprocessing

def get_cpu_count():
    return multiprocessing.cpu_count()

