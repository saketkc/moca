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
