import os
import subprocess
import time
from moca.helpers.encode import *

__root_dir__ = '/media/data1/encode_analysis'

records = get_idr_controlled_peaks()
for record in records:
    info = download_idr_tfs(__root_dir__, record)
    if info['assembly'] == 'hg19':
        cmd = 'mocacli -i {} -c application.cfg -g {} -o {}'.format(info['bedfile'], info['assembly'], os.path.join(os.path.dirname(info['bedfile']), 'moca_output'))
        start = time.time()
        proc = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        end = time.time()
        print 'stderr: {}'.format(stderr)
        print '{} sec'.format(end-start)
        if stderr:
            with open('errors.txt', 'a') as f:
                f.write(cmd)
        else:
            with open('processed.txt', 'a') as f:
                f.write(cmd)



