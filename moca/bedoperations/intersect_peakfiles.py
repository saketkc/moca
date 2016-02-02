import os
from format_peakfile import convert_to_scorefile
import ntpath
import sys
import pandas

broadPeak_columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'p-value', 'q-value']
narrowPeak_columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'p-value', 'q-value', 'peak']
questPeak_columns = ['chrom', 'position', 'score']

DIFFERENCE_THRESHOLD = 20
COUNT_THRESHOLD = 2
class Site:
    def __init__(self, chrom,  position, score):
        self.chrom = chrom
        self.position = position
        self.score = score
        self.lower_pos = self.position-DIFFERENCE_THRESHOLD
        self.upper_pos =  self.position+DIFFERENCE_THRESHOLD

def is_site_in_limit(source, target):
    return (target.position >= source.lower_pos) and (target.position <= source.upper_pos)

def is_inrange(query, target):
    lower = query-DIFFERENCE_THRESHOLD
    upper = query+DIFFERENCE_THRESHOLD
    return (target>=lower)  and (target<=upper)

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def get_sorted_filename(f):
    filename, extension = os.path.splitext(f)
    fileout = filename+'_sorted'+extension
    return fileout

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

def batch_convert_to_scorefile(list_of_files):
    for query_file in list_of_files:
        filetype = determine_filetype(query_file)
        ## Assume same filetype for all other files
        query_df = pandas.read_table(query_file)
        if filetype == 'broadPeak':
            query_df.columns = broadPeak_columns
        elif filetype == 'narrowPeak':
            query_df.columns = narrowPeak_columns
        elif filetype == 'questPeak':
            query_df.columns = questPeak_columns

        if filetype in ['broadPeak', 'narrowPeak']:
            outfile = get_sorted_filename(os.path.abspath(query_file))
            convert_to_scorefile(query_file, filetype, outfile)

def account_keeping():
    """
    For any site with a higher index, its sites which belong to a the 'positions'
    in a lower index, should automatically be added
    """
    pass

def exclude(list_d, index):
    return [v for i,v in enumerate(list_d) if i not in index]
def find_intersecting_peaks(list_of_files):
    all_sites_list_map = []
    all_chr = []
    for f in list_of_files:
        with open(f,'r') as fh:
            mapper = {}
            for line in fh:
                split = line.split('\t')
                chrom = split[0]
                position = int(split[1])
                score = float(split[2])
                if chrom not in mapper.keys():
                    mapper[chrom] = {}
                if chrom not in all_chr:
                    all_chr.append(chrom)
                assert position not in mapper[chrom].keys()
                mapper[chrom][position] = {'positions':[], 'indices':[], 'score' : score}
            all_sites_list_map.append(mapper)

    master = []
    for index, query_site_dict in enumerate(all_sites_list_map):
        pointer = index
        for target_site_dict in all_sites_list_map[index+1:]:
            pointer +=1
            for chrom in query_site_dict.keys():
                try:
                    target_chrom_positions = sorted(target_site_dict[chrom].keys())
                    for query_pos in query_site_dict[chrom].keys():
                        for target_pos in target_chrom_positions:
                            if is_inrange(query_pos, target_pos):
                                query_site_dict[chrom][query_pos]['positions'].append(target_pos)
                                query_site_dict[chrom][query_pos]['indices'].append(pointer)
                                all_sites_list_map[index] = query_site_dict
                except KeyError:
                    pass
        for j, previous in enumerate(master):
            for chrom in query_site_dict.keys():
                try:
                    target_chrom_positions = sorted(previous[chrom].keys())
                    for query_pos in query_site_dict[chrom].keys():
                        for target_pos in target_chrom_positions:
                            if is_inrange(query_pos, target_pos):
                                query_site_dict[chrom][query_pos]['positions'].append(target_pos)
                                query_site_dict[chrom][query_pos]['indices'].append(j)
                                all_sites_list_map[index] = query_site_dict
                except KeyError:
                    pass

        master.append(all_sites_list_map[index])

    for i,row in enumerate(master):
        print "######## Query {}############".format(i)
        for chrom in row.keys():
            for pos in sorted(row[chrom].keys()):
                print chrom, pos, row[chrom][pos]['positions'], row[chrom][pos]['indices']


    ## Start outputing based on threshold based cutoof
    out_file = 'common_intersection.bed'
    content = ''
    for i,row in enumerate(master):
        for chrom in row.keys():
            for pos in row[chrom].keys():
                indices = row[chrom][pos]['indices']
                target_pos = row[chrom][pos]['positions']
                scores = [master[indices[j]][chrom][target_pos[j]]['score'] for j in range(0, len(indices))]
                scores.append(row[chrom][pos]['score'])
                if len(scores) > COUNT_THRESHOLD:
                    score = sum(scores)/len(scores)*1.0
                    content+='{}\t{}\t{}\n'.format(chrom, pos, score)

    with open(out_file, 'w') as fp:
        fp.write(content)

    ##TODO Perform chromosome wise interesectipons




if __name__ == '__main__':
    list_of_files = sys.argv[1:]
    #batch_convert_to_scorefile(list_of_files)
    find_intersecting_peaks(list_of_files)
