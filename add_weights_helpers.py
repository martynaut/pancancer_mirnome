import requests
from xml.etree import ElementTree
import regex

rev_comp = {
    'G': 'C',
    'C': 'G',
    'A': 'U',
    'T': 'A',
    'U': 'A'
}

iupac_dict = {
    'R': '[GA]{1}',
    'Y': '[UC]{1}',
    'M': '[AC]{1}',
    'K': '[GU]{1}',
    'S': '[GC]{1}',
    'W': '[AU]{1}',
    'H': '[ACU]{1}',
    'B': '[GUC]{1}',
    'V': '[GCA]{1}',
    'D': '[GUA]{1}',
    'N': '[GCUA]{1}',
}

motifs = {
    'UAGGGAW': 'loop',
    'AUUUUUAUUUU': 'loop',
    'GGGU': 'loop',
    'GGAG': 'loop',
    'YGCY': 'loop',
    'UGC': 'loop',
    'UGU': 'loop',
    'AUCUU': 'loop',
    'SMUANY': 'loop',
    'CAUC': 'loop',
    'UAUAA': 'loop',
    'UUUUUCCNUCUUU': 'loop',
    'VCAUCH': 'all',
    'GCAUG': 'all',
    'CAGAC': 'all',
    'UGUNNNNNNNUGU': 'all',
    'GCAGCGC': 'all'
}


def mirna_weight_calc(row):
    all_nt = len(row['whole_seq'])
    if row['balance'] == '5p' or row['balance'] == '3p':
        nt_with_2 = 7
    else:
        nt_with_2 = 14
    nt_cut_only = 4

    if row['orientation'] == '+' and row['balance'] == '5p':
        nt_duplex = row[('stop', 'post-seed', '5p')] - row[('start', 'pre-seed', '5p')] + 1 - 7 +\
                    row[('stop', 'post-seed', '3p')] - row[('start', 'pre-seed', '3p')] + 1
    elif row['orientation'] == '+' and row['balance'] == '3p':
        nt_duplex = row[('stop', 'post-seed', '5p')] - row[('start', 'pre-seed', '5p')] + 1 +\
                    row[('stop', 'post-seed', '3p')] - row[('start', 'pre-seed', '3p')] + 1 - 7
    elif row['orientation'] == '+' and ((row['balance'] == 'both') or (row['balance'] == 'unknown')):
        nt_duplex = row[('stop', 'post-seed', '5p')] - row[('start', 'pre-seed', '5p')] + 1 - 7 +\
                    row[('stop', 'post-seed', '3p')] - row[('start', 'pre-seed', '3p')] + 1 - 7
    elif row['orientation'] == '-' and row['balance'] == '5p':
        nt_duplex = row[('stop', 'pre-seed', '5p')] - row[('start', 'post-seed', '5p')] + 1 - 7 +\
                    row[('stop', 'pre-seed', '3p')] - row[('start', 'post-seed', '3p')] + 1
    elif row['orientation'] == '-' and row['balance'] == '3p':
        nt_duplex = row[('stop', 'pre-seed', '5p')] - row[('start', 'post-seed', '5p')] + 1 +\
                    row[('stop', 'pre-seed', '3p')] - row[('start', 'post-seed', '3p')] + 1 - 7
    elif row['orientation'] == '-' and ((row['balance'] == 'both') or (row['balance'] == 'unknown')):
        nt_duplex = row[('stop', 'pre-seed', '5p')] - row[('start', 'post-seed', '5p')] + 1 - 7 +\
                    row[('stop', 'pre-seed', '3p')] - row[('start', 'post-seed', '3p')] + 1 - 7
    all_motifs = []
    for motif in motifs.keys():
        regex_motif = motif
        for letter, reg in iupac_dict.items():
            regex_motif = regex_motif.replace(letter, reg)
        all_before = regex.finditer(regex_motif, row['whole_seq'], overlapped=True)
        all_loc_loop = []
        all_loc_flanks = []
        for i in all_before:

            motif_start, motif_stop = i.span()

            if row['orientation'] == '+':
                all_loc = [x + row[1] + 1 for x in range(motif_start, motif_stop, 1)]
                all_loc_loop = [x for x in all_loc if (x > row[('stop', 'post-seed', '5p')] + 1) and
                                (x < row[('start', 'pre-seed', '3p')] - 1)]
                all_loc_flanks = [x for x in all_loc if (x < row[('start', 'pre-seed', '5p')] - 1) or
                                (x > row[('stop', 'post-seed', '3p')] + 1)]
            elif row['orientation'] == '-':
                all_loc = [row[2] - x for x in range(motif_start, motif_stop, 1)]
                all_loc_loop = [x for x in all_loc if (x < row[('start', 'post-seed', '5p')] - 1) and
                                (x > row[('stop', 'pre-seed', '3p')] + 1)]
                all_loc_flanks = [x for x in all_loc if (x > row[('stop', 'pre-seed', '5p')] + 1) or
                                  (x < row[('start', 'post-seed', '3p')] - 1)]

            if motifs[motif] == 'loop':
                all_motifs = all_motifs + all_loc_loop
            else:
                all_motifs = all_motifs + all_loc_loop + all_loc_flanks
    all_motifs = sorted(list(set(all_motifs)))
    len_motifs = len(all_motifs)
    weight_cal = all_nt - (nt_with_2+len_motifs+nt_duplex+nt_cut_only) + nt_with_2*2 + \
                 (len_motifs+nt_duplex+nt_cut_only)* 1.5

    return weight_cal


def motif_check(row_values, regex_motif, raw_motif):
    if motifs[raw_motif] == 'loop':
        if row_values['type'] != 'loop':
            return 'no_change'
    all_before = regex.finditer(regex_motif, row_values['whole_seq'], overlapped=True)
    for i in all_before:
        motif_start, motif_stop = i.span()
        if row_values['orientation'] == '+':
            if (row_values['pos'] < motif_stop + row_values[1]) \
                    and (row_values['pos'] > motif_start + row_values[1]):
                return 'loss'
        else:
            if (row_values['pos'] < row_values[2] - motif_start) \
                    and (row_values['pos'] > row_values[2] - motif_stop):
                return 'loss'
    return 'no_change'


def change_in_motifs(row_value, motifs_columns):
    for col in motifs_columns:
        if row_value[col] != 'no_change':
            return 1
    return 0


def get_whole_sequence(row):
    URL = "http://genome.ucsc.edu/cgi-bin/das/hg38/dna"
    PARAMS = {'segment': str(row[0]) + ':' + str(row[1]+1) + ',' + str(row[2])}
    r = requests.get(url=URL, params=PARAMS)
    tree = ElementTree.fromstring(r.content)
    seq = []
    for child in tree:
        for cc in child:
            seq.append(cc.text.upper())
    if len(seq) == 1:
        return seq[0].replace(' ', '').replace('\n', '').strip()
    else:
        raise ValueError


def change_sequence(row):
    seq = row['whole_seq']
    if ',' in row['alt']:
        row['alt'] = row['alt'].split(',')[0]
    if row['orientation'] == '+':
        pos = row['pos'] - row[1]
        ref = row['ref'].replace('T', 'U')
        alt = row['alt'].replace('T', 'U')
        if len(ref) > len(seq[pos:pos + len(ref)]) and \
                (pos + len(ref)) > 0:
            ref = ref[:len(seq[pos:pos + len(ref)])]
        if seq[pos:pos + len(ref)] == ref:
            seq_mut = seq[:pos] + alt + seq[pos + len(ref):]
        else:
            raise ValueError
    else:
        pos = row[2] - row['pos']
        ref = ''.join([rev_comp[c] for c in list(row['ref'])])[::-1]
        alt = ''.join([rev_comp[c] for c in list(row['alt'])])[::-1]
        if len(ref) > len(seq[max([pos - len(ref) +1, 0]):pos+1]) and \
                (pos - len(ref) + 1) < 0:
            ref = ref[len(ref) - len(seq[max([pos - len(ref) +1, 0]):pos+1]):]
        elif len(ref) > len(seq[max([pos - len(ref) +1, 0]):pos+1]) and \
                (pos - len(ref) + 1) >= 0:
            ref = ref[:len(seq[pos - len(ref) +1:pos+1])]
        if seq[pos - len(ref) +1:pos+1] == ref:
            seq_mut = seq[:pos] + alt + seq[pos + len(ref):]
        else:
            raise ValueError

    return seq_mut
