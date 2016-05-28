"""Utility function to run jobs
"""
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import contextlib
import os
import re
import subprocess
import time

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.

    Credits: Brad Chapman for bcbio-nextgen: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/utils.py#L172

    Parameters
    ----------
    dname: str
        Path of directory to be created
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    Credits: Brad Chapman for bcbio-nextgen: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/utils.py#L192

    Parameters
    ----------
    new_dir: str
        Location of directory to be created
    """
    cur_dir = os.getcwd()
    safe_makedir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)

def run_job(cmd, cwd):
    """Execute command line jobs passed as str argument

    Parameters
    ----------
    cmd: str
        Absolute command line statement as would be run on command line

    Returns
    -------
    stdout: str
        standard output as written by cmd program
    stderr: str
        standard output as written by cmd program
    returncode: int
        return code as returned by subprocess(not very informative in general)
    """

    # Replace multiple spaces with single space
    cmd = re.sub('\s+', ' ', cmd).strip()
    cmd_split = cmd.split(' ')
    with chdir(cwd):
        proc = subprocess.Popen(cmd_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    returncode = proc.returncode
    return stdout, stderr, returncode
