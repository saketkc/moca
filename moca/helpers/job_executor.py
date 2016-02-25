"""Utility function to run jobs
"""

import subprocess
import os
import contextlib
import time


def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    Credits: Brad Chapman for bcbio-nextgen:
        https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/utils.py#L172
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

    Credits: Brad Chapman for bcbio-nextgen:
        https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/utils.py#L192
    """
    cur_dir = os.getcwd()
    safe_makedir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)

def run_job(cmd, cwd):
    """Execute command line jobs passed as string argument

    Parameters
    ----------
    cmd: string
        Absolute command line statement as would be run on command line

    Returns
    -------
    stdout: string
        standard output as written by cmd program
    stderr: string
        standard output as written by cmd program
    returncode: int
        return code as returned by subprocess(not very informative in general)
    """

    cmd_split = cmd.split(' ')
    with chdir(cwd):
        proc = subprocess.Popen(cmd_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    returncode = proc.returncode
    return stdout, stderr, returncode
