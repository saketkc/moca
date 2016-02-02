#!/usr/bin/env python
import pandas
import sys
import os
import shutil

column_names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'p-value', 'q-value']

def determine_filetype(filepath):
    p = pandas.read_table(filepath)
    c = len(p.columns)
    if c==10:
        return 'narrowPeak'
    elif c==9:
        return 'broadPeak'
    elif c==3:
        return 'questPeak'
    return 'unknown'

def convert_to_scorefile(filepath, filetype=None, out_file=None):
    if not filetype:
        filetype = determine_filetype(filepath)
    assert filetype in ['broadPeak', 'narrowPeak', 'questPeak']
    if filetype == 'questPeak':
        ## simply copy to peaks.tsvi
        assert out_file is not None
        if filepath!=out_file:
            # Copy only ifthey do not have samename
            shutil.copy(filepath, out_file)
        return True

    df = pandas.read_table(filepath, header=None)
    ## if file is narrow peak, the speak summit occurs
    #3 at 0-based offset given in the last columi
    if filetype=='broadPeak':
        df.columns = column_names
    else:
        columns = column_names[:]
        columns.append('peak')
        df.columns = columns

    if filetype=='narrowPeak':
        filter_df1 = df[df.peak.astype(int)==-1]
        filter_df2 = df[df.peak.astype(int)!=-1]
        filter_df1['peak_positions'] = (filter_df1['chromStart'].astype(int)+filter_df1['chromEnd'].astype(int))
        filter_df1['peak_positions'] = [int(x/2) for x in filter_df1['peak_positions'].astype(int)]
        filter_df2['peak_positions'] = filter_df2['chromStart'].astype(int)+filter_df2['peak'].astype(int)
        df = pandas.concat([filter_df1, filter_df2])
    else:
        df['peak_positions'] = (df['chromStart']+df['chromEnd'])
        df['peak_positions'] = [int(x/2) for x in df['peak_positions'].astype(int)]

    if not out_file:
        out_file = os.path.join(os.path.dirname(filepath), 'peaks.tsv')
    df = df.sort(columns=['score'], ascending=False)
    df['peak_positions'] = df['peak_positions'].astype(int)
    df.to_csv(out_file, sep='\t', columns=['chrom', 'peak_positions', 'score'], index=False, header=False)

if __name__ == '__main__':
    convert_to_scorefile(sys.argv[1], sys.argv[2], sys.argv[3])

