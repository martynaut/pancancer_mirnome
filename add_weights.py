import pandas as pd
import scipy.stats as ss
from add_weights_helpers import rev_comp, change_in_motifs, iupac_dict, motifs, \
    motif_check, mirna_weight_calc
# from add_weights_helpers import get_whole_sequence


df_all_mut = pd.read_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/all_mutations.csv'
)

df_all_mut.drop('Unnamed: 0', axis=1, inplace=True)

print(df_all_mut.head())

df_patient_count = pd.read_excel(
    '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/patient_mutation_count.xlsx'
)

df_all_mut = df_all_mut.join(df_patient_count.set_index(['patient_id', 'cancer_id']),
                             ['indiv_name', 'cancer'],
                             how='left')

df_all_mut['>10000'] = df_all_mut['mutation_count'].apply(
    lambda x: 1 if x > 10000 else 0
)

df_all_mut['cut_region'] = df_all_mut.apply(
    lambda x: 1 if
    ((x['type'] == 'flanking-5' and x['from end'] == -1)
        or (x['type'] == 'flanking-3' and x['from_start'] == 1)
        or (x['type'] == 'loop' and x['from_start'] == 1)
        or (x['type'] == 'loop' and x['from end'] == -1)
        or (x['type'] == 'pre-seed' and x['from_start'] == 1)
        or (x['type'] == 'post-seed' and x['from end'] == -1)
     ) else 0,
    axis=1
)

df_all_mut['duplex'] = df_all_mut['type'].apply(
    lambda x: 1 if x in ['pre-seed', 'seed', 'post-seed'] else 0
)

df_all_mut['seed'] = df_all_mut.apply(
    lambda x: 1 if x['type'] == 'seed' and ((x['arm'] == x['balance'])
                                            or x['balance'] == 'both') else 0,
    axis=1
)

# df_coordinates = pd.read_csv(
#     '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates.bed',
#     sep='\t',
#     header=None
# )
#
# df_coordinates['whole_seq'] = df_coordinates.apply(
#     get_whole_sequence,
#     axis=1
# )
#
# df_coordinates.to_csv(
#     '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates_withseq.bed',
#     sep='\t',
#     index=False,
#     header=False
# )

df_coordinates = pd.read_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates_withseq.bed',
    sep='\t',
    header=None
)

print(df_coordinates.head())

print(df_coordinates[df_coordinates.duplicated(subset=[0, 3],
                                               keep='first')])

print(df_all_mut.shape)

df_all_mut = df_all_mut.join(df_coordinates.set_index([0, 3]),
                             on=['chrom', 'pre_name'],
                             how='left')

df_all_mut = df_all_mut[(df_all_mut['pos'] >= df_all_mut[1])
                        & (df_all_mut['pos'] <= df_all_mut[2])]

print(df_all_mut.shape)

df_all_mut['whole_seq'] = df_all_mut[4].str.upper()

df_all_mut['whole_seq'] = df_all_mut.apply(
    lambda x: x['whole_seq'].replace('T', 'U') if x['orientation'] == '+'
    else ''.join([rev_comp[c] for c in list(x['whole_seq'])])[::-1],
    axis=1
)

# df_all_mut['mutation_seq'] = df_all_mut.apply(
#     change_sequence,
#     axis=1
# )


for key in motifs.keys():
    motif = key
    for letter, reg in iupac_dict.items():
        motif = motif.replace(letter, reg)
    df_all_mut['motifs_{}'.format(key)] = df_all_mut.apply(lambda x: motif_check(x, motif, key), axis=1)

motifs_cols = [col for col in df_all_mut.columns if 'motifs_' in str(col)]
print(motifs_cols)

df_all_mut.head()

df_all_mut['motifs'] = df_all_mut.apply(
    lambda x: change_in_motifs(x, motifs_cols),
    axis=1
)

df_all_mut['weight'] = df_all_mut.apply(
    lambda x: 2 if x['seed'] == 1 else (
        1.5 if (x['duplex'] == 1)
        or (x['motifs'] == 1) or (x['cut_region'] == 1) else 1
    ),
    axis=1
)

print(df_all_mut.head())

df_all_mut.to_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/all_mutations_weights.csv'
)

df_localization = pd.read_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/localizations.csv'
)

df_localization = df_localization.groupby(['chrom', 'pre_name', 'start_pre', 'orientation',
                                           'balance', 'type', 'arm']).min()
df_localization = df_localization.unstack(5).unstack(5)[['start', 'stop']]
print(df_localization.head())
print(df_localization.shape)
df_localization = df_localization.join(df_coordinates.set_index([0, 3]), ['chrom', 'pre_name'], 'inner')
print(df_localization.shape)
df_localization = df_localization.reset_index()
df_localization = df_localization[(df_localization[('start', 'seed', '5p')] <= df_localization[2])
                                  & (df_localization[('start', 'seed', '5p')] >= df_localization[1])]
print(df_localization.head())
print(df_localization.shape)
df_localization['whole_seq'] = df_localization[4].str.upper()

df_localization['whole_seq'] = df_localization.apply(
    lambda x: x['whole_seq'].replace('T', 'U') if x['orientation'] == '+'
    else ''.join([rev_comp[c] for c in list(x['whole_seq'])])[::-1],
    axis=1
)
df_localization['weight'] = df_localization.apply(lambda x: mirna_weight_calc(x), axis=1)

df_localization.to_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/localizations_weight.csv'
)

df_localization['weight_prob'] = df_localization['weight'] / df_localization['weight'].sum()


coords = pd.read_csv(
            '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates.bed',
            sep='\t', header=None)
coords['length'] = coords[2] - coords[1]
all_mirnas_len = coords[coords[3].str.contains('hsa-')]['length'].sum()
coords['prob'] = coords['length'] / all_mirnas_len


df_occur = pd.DataFrame()
for cancer in df_all_mut['cancer'].unique():
    df_temp_raw = df_all_mut[df_all_mut['cancer'] == cancer]

    no_patients_before = df_temp_raw['indiv_name'].nunique()

    df_temp_raw = df_temp_raw[df_temp_raw['>10000'] != 1]

    no_patients = df_temp_raw['indiv_name'].nunique()

    no_of_mutations = df_temp_raw.shape[0]
    no_of_pos = df_temp_raw.groupby(['indiv_name', 'chrom', 'pos']).count().shape[0]

    if df_temp_raw.shape[0] > 0:

        df_temp = df_temp_raw[['pre_name', 1, 'indiv_name']].groupby(['pre_name', 1], as_index=False).count()

        mutation_count = df_temp['indiv_name'].sum()

        df_temp['binom'] = df_temp.apply(lambda x: ss.binom_test(x['indiv_name'],
                                                                 mutation_count,
                                                                 coords[coords[3] == x['pre_name']]['prob'].values[0],
                                                                 alternative='greater'), axis=1)

        df_temp['order'] = df_temp['binom'].apply(lambda x: df_temp[df_temp['binom'] <= x].shape[0]
                                                  )
        df_temp.sort_values('order', inplace=True)
        df_temp['FDR'] = df_temp.apply(lambda x: x['order'] / 1918 * 0.01, axis=1)
        df_temp['signif'] = df_temp['binom'] < df_temp['FDR']
        # print(cancer)
        # print(df_temp[['pre_name', 'indiv_name', 'binom', 'order', 'FDR', 'signif']].head(20))
        max_order = df_temp[df_temp['signif']]['order'].max()
        # print(max_order)
        df_temp['BH'] = df_temp.apply(lambda x: x['binom'] * 1918 / x['order'], axis=1)
        mirnas = df_temp.loc[df_temp['order'] <= max_order, 'pre_name'].unique()

        df_temp2 = df_temp_raw[['pre_name', 1, 'weight']].groupby(['pre_name', 1], as_index=False).sum()

        mutation_count2 = df_temp2['weight'].sum()

        df_temp2['binom'] = df_temp2.apply(lambda x: ss.binom_test(x['weight'],
                                                                   mutation_count2,
                                                                   df_localization[
                                                                       df_localization['pre_name'] == x['pre_name']][
                                                                           'weight_prob'].values[0],
                                                                   alternative='greater'), axis=1)

        df_temp2['order'] = df_temp2['binom'].apply(lambda x: df_temp2[df_temp2['binom'] <= x].shape[0])
        df_temp2.sort_values('order', inplace=True)
        df_temp2['FDR'] = df_temp2.apply(lambda x: x['order'] / 1918 * 0.01, axis=1)
        df_temp2['signif'] = df_temp2['binom'] < df_temp2['FDR']
        max_order2 = df_temp2.loc[df_temp2['signif'], 'order'].max()
        df_temp2['BH'] = df_temp2.apply(lambda x: x['binom'] * 1918 / x['order'], axis=1)
        mirnas_weight = df_temp2.loc[df_temp2['order'] <= max_order, 'pre_name'].unique()

        mutated_gene_count = df_temp.shape[0]

    else:
        mirnas = []
        mirnas_weight = []
        mutated_gene_count = 0

    mirnas_count = len(mirnas)
    mirnas_weight_count = len(mirnas_weight)
    row = pd.DataFrame([[cancer,
                         no_patients_before,
                         no_patients,
                         mutated_gene_count,
                         no_of_mutations,
                         no_of_pos,
                         mirnas, mirnas_count,
                         mirnas_weight, mirnas_weight_count
                         ]], columns=['cancer',
                                      'no_patients_before',
                                      'no_patients',
                                      'mutated_genes_count',
                                      'no_of_mutations',
                     'no_of_pos',
                                      'significant_mirnas',
                                      'significant_mirnas_count',
                                      'significant_mirnas_weighted',
                                      'significant_mirnas_weighted_count'])
    df_occur = pd.concat([df_occur, row])

df_occur = df_occur.sort_values('cancer', ascending=True)

no_patients_before = df_all_mut['indiv_name'].nunique()
df_temp_raw = df_all_mut[df_all_mut['>10000'] != 1]
no_patients = df_temp_raw['indiv_name'].nunique()

no_of_mutations = df_temp_raw.shape[0]
no_of_pos = df_temp_raw.groupby(['indiv_name', 'chrom', 'pos']).count().shape[0]

if df_temp_raw.shape[0] > 0:

    df_temp = df_temp_raw[['pre_name', 1, 'indiv_name']].groupby(['pre_name', 1], as_index=False).count()

    mutation_count = df_temp['indiv_name'].sum()

    df_temp['binom'] = df_temp.apply(lambda x: ss.binom_test(x['indiv_name'],
                                                             mutation_count,
                                                             coords[coords[3] == x['pre_name']]['prob'].values[0],
                                                             alternative='greater'), axis=1)

    df_temp['order'] = df_temp['binom'].apply(lambda x: df_temp[df_temp['binom'] <= x].shape[0]
                                              )
    df_temp.sort_values('order', inplace=True)
    df_temp['FDR'] = df_temp.apply(lambda x: x['order'] / 1918 * 0.01, axis=1)
    df_temp['signif'] = df_temp['binom'] < df_temp['FDR']
    # print(cancer)
    # print(df_temp[['pre_name', 'indiv_name', 'binom', 'order', 'FDR', 'signif']].head(20))
    max_order = df_temp[df_temp['signif']]['order'].max()
    # print(max_order)
    df_temp['BH'] = df_temp.apply(lambda x: x['binom'] * 1918 / x['order'], axis=1)
    mirnas = df_temp.loc[df_temp['order'] <= max_order, 'pre_name'].unique()

    df_temp2 = df_temp_raw[['pre_name', 1, 'weight']].groupby(['pre_name', 1], as_index=False).sum()

    mutation_count2 = df_temp2['weight'].sum()

    df_temp2['binom'] = df_temp2.apply(lambda x: ss.binom_test(x['weight'],
                                                               mutation_count2,
                                                               df_localization[
                                                                   df_localization['pre_name'] == x['pre_name']][
                                                                       'weight_prob'].values[0],
                                                               alternative='greater'), axis=1)

    df_temp2['order'] = df_temp2['binom'].apply(lambda x: df_temp2[df_temp2['binom'] <= x].shape[0])
    df_temp2.sort_values('order', inplace=True)
    df_temp2['FDR'] = df_temp2.apply(lambda x: x['order'] / 1918 * 0.01, axis=1)
    df_temp2['signif'] = df_temp2['binom'] < df_temp2['FDR']
    max_order2 = df_temp2.loc[df_temp2['signif'], 'order'].max()
    df_temp2['BH'] = df_temp2.apply(lambda x: x['binom'] * 1918 / x['order'], axis=1)
    mirnas_weight = df_temp2.loc[df_temp2['order'] <= max_order, 'pre_name'].unique()

    mutated_gene_count = df_temp.shape[0]

else:
    mirnas = []
    mirnas_weight = []
    mutated_gene_count = 0

mirnas_count = len(mirnas)
mirnas_weight_count = len(mirnas_weight)
row = pd.DataFrame([['pancancer',
                     no_patients_before,
                     no_patients, mutated_gene_count,
                     no_of_mutations,
                     no_of_pos,
                     mirnas, mirnas_count,
                     mirnas_weight, mirnas_weight_count
                     ]], columns=['cancer',
                                  'no_patients_before',
                                  'no_patients',
                                  'mutated_genes_count',
'no_of_mutations',
                     'no_of_pos',
                                  'significant_mirnas',
                                  'significant_mirnas_count',
                                  'significant_mirnas_weighted',
                                  'significant_mirnas_weighted_count'])
df_occur = pd.concat([df_occur, row])

df_occur.set_index('cancer', inplace=True)
df_occur.to_excel('/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/weight_summary.xlsx')
