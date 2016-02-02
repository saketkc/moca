#!/usr/bin/env python
"""
Script to check if sites outputted by fimo
are overlapping.

Overlap is defined by this example
chr1 12 20
chr1 13 21

Overlap: 13-20
"""
import pandas
import sys
__column_names__ = ['pattern_name','sequence_name','start','stop','strand','score','p_value', 'q_value', 'matched_sequence']
def fimo_sites_intersect(fimo_file):
    table = pandas.read_table(fimo_file)
    table.columns = __column_names__
    coord_fields = pandas.DataFrame(table['sequence_name'].str.split('_').tolist(), columns = ['chr', 'start', 'end'])
    table['chr'] = coord_fields['chr']
    table['chromStart'] = coord_fields['start'].astype(int)+coord_fields['end'].astype(int)+table['start'].astype(int)
    table['chromEnd'] = coord_fields['start'].astype(int)+coord_fields['end'].astype(int)+table['stop'].astype(int)+1
    table = table.sort(['chr'])
    table.to_csv('fimo_2_sites.txt', sep='\t', index=False, columns=['chr', 'chromStart', 'chromEnd', 'strand', 'p_value', 'q_value'])

    fimo_sites = table[['chr', 'chromStart', 'chromEnd', 'strand', 'p_value', 'q_value', 'matched_sequence']]
    fimo_sites_sorted = fimo_sites.groupby('chr')
    for chrom,rows in fimo_sites_sorted.groups.iteritems():
        for index, outer_row in enumerate(rows):
            outer_row_content = fimo_sites.loc[[outer_row]]
            outer_row_series = set(range(outer_row_content['chromStart'], outer_row_content['chromEnd']+1))
            outer_strand = outer_row_content['strand'].values[0]
            outer_match_seq = outer_row_content['matched_sequence'].values[0]
            for inner_row in rows[index+1:]:
                inner_row_content = fimo_sites.loc[[inner_row]]
                inner_strand = inner_row_content['strand'].values[0]
                inner_row_series = set(range(inner_row_content['chromStart'], inner_row_content['chromEnd']+1))
                inner_match_seq = inner_row_content['matched_sequence'].values[0]
                assert(len(inner_row_series)==len(outer_row_series))
                intersection = outer_row_series.intersection(inner_row_series)
                if intersection and (outer_strand == inner_strand):
                    print '{}: {} {} {} {} {} {}\n'.format(chrom,
                                                           outer_row_content['chromStart'].values[0],
                                                           outer_row_content['chromEnd'].values[0],
                                                           inner_row_content['chromStart'].values[0],
                                                           inner_row_content['chromEnd'].values[0],
                                                           outer_match_seq,
                                                           inner_match_seq)

if __name__ == '__main__':
    fimo_sites_intersect(sys.argv[1])
