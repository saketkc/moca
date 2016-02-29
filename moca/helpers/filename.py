import ntpath
import os

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def filename_extension(path):
    fn = path_leaf(path)
    return os.path.splitext(fn)

