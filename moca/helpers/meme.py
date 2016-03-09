from Bio import motifs

def read_memefile(meme_file):
    """
    Read meme file
    """
    summary = {}
    summary['motif_occurrences'] = {}
    records = motifs.parse(open(meme_file), 'meme')
    summary['total_motifs'] = len(records)
    num_occurrences = []
    for index, record in enumerate(records):
        num_occurrences.append(int(getattr(record,'num_occurrences','Unknown')))

    sorted_occurences = sorted(enumerate(num_occurrences), key=lambda x: x[1])
    summary['motif_occurrences'] = {'motif{}'.format(index+1):value for index,value in sorted_occurences}
    summary['motif_records'] = records
    return summary

def position_wise_profile(counts_dict, length):
    """
    Convert base to position wise profile
    """
    profile = map(dict, zip(*[[(k, v) for v in value] for k, value in counts_dict.items()]))
    return profile

def find_max_occurence(profile, max_count=2):
    """
    Return profile with base corresponding to max scores[CHECK!]
    """
    sorted_profile = []
    for p in profile:
        sorted_profile.append(sorted(p.items(), key=lambda x:x[1]))
    for i,p in enumerate(sorted_profile):
        sorted_profile[i] = p[-max_count:]
    return sorted_profile
