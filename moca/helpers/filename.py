from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
import ntpath
import fnmatch
import os

def touch(fp):
    """Touch an empty file"""
    open(fp, 'a').close()

def path_leaf(path):
    """Return head, tail of an absolute path"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def filename_extension(path):
    """Return filename, extension tuple"""
    fn = path_leaf(path)
    return os.path.splitext(fn)

def get_filename_without_ext(path):
    """Get filename from full path
    Parameters
    ---------
    path: str
        Absolute path

    """
    return os.path.splitext(path_leaf(path))[0]

def search_files(root_dir, ext):
    """Search files recursively with certain extension
    Credits: http://stackoverflow.com/a/2186565/756986
    """
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, '*{}'.format(ext) ):
            matches.append(os.path.join(root, filename))
    return matches

