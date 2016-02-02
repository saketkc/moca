import sys
import os
import subprocess
import ntpath
import pandas
from format_peakfile import convert_to_scorefile
import numpy as np

OVERLAP_FRACTION = 0.1
column_names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'p-value', 'q-value']

def is_broad_or_narrow(filename):
    p = pandas.read_table(filename)
    c = len(p.columns)
    if c==10:
        return 'narrowPeak'
    elif c==9:
        return 'broadPeak'
    return 'unknown'

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
def run_subprocess(command, cwd=None, stdout=None):
    print command
    split_command = command.split(' ')
    output = subprocess.call(split_command, cwd=cwd, stdout=stdout)
    if (int(output)!=0):
        sys.stderr.write('Error executing: {}'.format(command))
        sys.exit(1)
    return output

def count_lines(file):
    with open(file) as f:
        for i, line in enumerate(f):
            pass
    return i+1

def sort_files_by_size(list_of_files):
    """
    Sort files in decending order of file size
    """
    filesizes = [count_lines(f) for f in list_of_files]
    sorted_filesizes = sorted(enumerate(filesizes), key=lambda x: x[1], reverse=True)
    sorted_filenames = [list_of_files[i] for i,v in sorted_filesizes]
    return sorted_filenames

def get_sorted_filename(f):
    filename, extension = os.path.splitext(f)
    fileout = filename+'_sorted'+extension
    return fileout

def sort_bedfiles(list_of_files):
    for f in list_of_files:
        fileout = get_sorted_filename(f)
        cwd = os.path.dirname(f)
        f_abs = path_leaf(f)
        handle = open(fileout,'w')
        run_subprocess('bedtools sort -i {}'.format(f_abs), stdout=handle, cwd=cwd)
        handle.close()

def filter_peaks(filepath, filetype, out_file):
    df = pandas.read_table(filepath, header=None)
    if filetype=='broadPeak':
        columns = column_names[:]
        columns.append('intersection_count')
        df.columns = columns
    else:
        columns = column_names[:]
        columns.append('peak')
        columns.append('intersection_count')
        df.columns = columns
    filter_df = df[df.intersection_count>=1]
    filter_df.to_csv(out_file, sep='\t', columns=columns[:-1], index=False, header=False)

def perform_bed_intersection(query_file, list_of_target_files):
    target_files = [os.path.abspath(f) for f in list_of_target_files]
    target_files_abs = (' ').join(target_files)
    query_file_abs = os.path.abspath(query_file)
    command = 'bedtools intersect -c -sorted -sortout -a {} -b {}'.format(query_file_abs, target_files_abs)
    outfile = os.path.join(os.path.dirname(query_file), 'intersection.bed')
    handle = open(outfile, 'w')
    cwd = os.path.dirname(query_file)
    run_subprocess(command, stdout=handle, cwd=cwd)
    handle.close()
    broad_or_narrow = is_broad_or_narrow(query_file_abs)
    print broad_or_narrow
    filter_outfile = os.path.join(os.path.dirname(query_file), 'intersection_filtered.bed')
    filter_peaks(outfile, broad_or_narrow, filter_outfile)
    convert_to_scorefile(filter_outfile, broad_or_narrow)

if __name__ == '__main__':
    files = sys.argv[1:]
    sorted_files = sort_files_by_size(files)
    sort_bedfiles(sorted_files)
    sorted_bedfiles = [get_sorted_filename(f) for f in sorted_files]
    perform_bed_intersection(sorted_bedfiles[0], sorted_bedfiles[1:])

