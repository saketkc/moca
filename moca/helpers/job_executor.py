import subprocess
import os
from contextlib import contextmanager

@contextmanager

def run_job(self, cmd, cwd):
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
    proc = subprocess.Popen(cmd_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderrr = proc.communicate()
    returncode = proc.returncode
    return stdout, stderr, returncode


