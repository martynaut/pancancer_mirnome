import pandas as pd
from distinct_occure_helpers import concat_ints, type_of_mutation, concat_alg, seq_type, if_complex, \
    subst_type, from_end, from_start


def extract_mirnas_maf(input_file, output_folder, hg19_file):

    hg19_coordinates = pd.read_csv(hg19_file, sep='\t')
    print(hg19_coordinates.head())

    n = 0

    output_file_name = output_folder + '/' + 'temp_mirnas_' + input_file.split('/')[-1]
    with open(input_file, 'r+') as input_maf:
        with open(output_file_name, 'w+') as output_maf:
            for line in input_maf:
                if line.startswith('Hugo'):
                    output_maf.write(line)
                else:
                    splitted = line.split('\t')
                    chrom = splitted[1]
                    start = int(splitted[2])
                    stop = int(splitted[3])
                    if hg19_coordinates[((hg19_coordinates['CHROMOSOME'] == chrom) &
                                         (hg19_coordinates['START'] <= start) &
                                         (hg19_coordinates['END'] >= start))
                                        | ((hg19_coordinates['CHROMOSOME'] == chrom) &
                                        (hg19_coordinates['START'] <= stop) &
                                         (hg19_coordinates['END'] >= stop)
                                         )].shape[0] > 0:
                        output_maf.write(line)
                        n = n + 1
                        print(n)


def liftover_mutations(input_file, output_folder, hg19_file, hg38_file):
    output_file_name = output_folder + '/' + 'temp_mirnas_' + input_file.split('/')[-1]

    hg19_coordinates = pd.read_csv(hg19_file, sep='\t')
    print(hg19_coordinates.head())

    hg38_coordinates = pd.read_csv(hg38_file, sep='\t')
    print(hg38_coordinates.head())
    print(hg19_coordinates.shape)
    hg19_coordinates = hg19_coordinates.join(hg38_coordinates.set_index('ELEMENT'),
                                             on='ELEMENT', how='inner', rsuffix='_hg38')
    print(hg19_coordinates.shape)
    hg19_coordinates.rename(index={'hsa-mir-4477b-1': 'hsa-mir-4477b',
                       'hsa-mir-4477b-2': 'hsa-mir-4477b',
                       'hsa-mir-4477a-1': 'hsa-mir-4477b',
                       'hsa-mir-4477a-2': 'hsa-mir-4477b',
                       'hsa-mir-10401-1': 'hsa-mir-10401',
                       'hsa-mir-10401-2': 'hsa-mir-10401',
                       'hsa-mir-10401-3': 'hsa-mir-10401',
                       'hsa-mir-10401-4': 'hsa-mir-10401'}, inplace=True)

    mutations = pd.read_csv(output_file_name, sep='\t')
    print(mutations.head())
    print(mutations.shape)

    mutations = mutations[['Chromosome', 'Start_position', 'End_position', 'Strand',
                           'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                           'Tumor_Sample_Barcode', 'Project_Code']]

    mutations['Alt'] = mutations.apply(lambda x: x['Tumor_Seq_Allele1'] if x['Tumor_Seq_Allele1'] != x['Reference_Allele'] else
                                       x['Tumor_Seq_Allele2'] if x['Tumor_Seq_Allele2'] != x['Reference_Allele'] else
                                       'NA', axis=1)
    print(mutations.head())
    mutations = mutations.join(hg19_coordinates.set_index('CHROMOSOME'), on='Chromosome', how='inner')
    mutations = mutations[mutations['Start_position'] <= mutations['END']]
    mutations = mutations[mutations['Start_position'] >= mutations['START']]

    print(mutations.shape)
    print(mutations[mutations.duplicated(subset=['Chromosome', 'Start_position', 'End_position', 'Strand',
                           'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                           'Tumor_Sample_Barcode', 'Project_Code'], keep=False)])

    mutations['Start_position'] = mutations['Start_position'] - mutations['START'] + mutations['START_hg38']
    mutations['End_position'] = mutations['End_position'] - mutations['START'] + mutations['START_hg38']

    mutations.to_csv(output_folder + '/' + 'all_mutations.csv', index=False)


def dist_occur(output_folder, localization_file, coordinates_file):
    mutations = pd.read_csv(output_folder + '/' + 'all_mutations.csv')
    print(mutations.head())
    print(mutations.shape)
    print(mutations.columns)
    mutations['Chromosome'] = mutations['Chromosome'].apply(lambda x: 'chr' + x)
    mutations.drop(['End_position', 'Strand', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                    'START', 'ELEMENT', 'START_hg38', 'END_hg38',
                    'END', 'CHROMOSOME_hg38'
                    ], axis=1, inplace=True)
    mutations.rename({'Chromosome': 'chrom',
                      'Start_position': 'pos',
                      'Reference_Allele': 'ref',
                      'Alt': 'alt',
                      'Tumor_Sample_Barcode': 'indiv_name',
                      'Project_Code': 'cancer'},
                     axis='columns', inplace=True)

    mutations['mutation_type'] = mutations.apply(lambda x: type_of_mutation(x), axis=1)

    mutations.fillna(-1, inplace=True)

    all_mutations = mutations.drop_duplicates(subset=['chrom', 'pos',
                                                      'indiv_name', 'ref', 'alt',
                                                      'mutation_type', 'cancer'],
                                              keep='first')
    print(all_mutations.shape)
    all_mutations['type_of_subst'] = all_mutations.apply(lambda x: subst_type(x), axis=1)

    localizations = pd.read_csv(localization_file, sep=',')
    print(all_mutations.shape)
    all_mutations = all_mutations.join(localizations.set_index('chrom'), on='chrom', how='left')
    all_mutations = all_mutations[(all_mutations['pos'] >= all_mutations['start'])
                                  & (all_mutations['pos'] <= all_mutations['stop'])]

    print(all_mutations.shape)
    double_check = all_mutations.groupby(['chrom', 'pos', 'indiv_name'])[['ref']].count()
    double_check.columns = ['no_of_loc']
    all_mutations = all_mutations.join(double_check, on=['chrom', 'pos', 'indiv_name'], how='left')

    print(all_mutations.shape)

    all_mutations['seq_type'] = all_mutations['name'].apply(lambda x: seq_type(x))

    all_mutations['from_start'] = all_mutations.apply(lambda x: from_start(x, 'start', 'stop'), axis=1)
    all_mutations['from end'] = all_mutations.apply(lambda x: from_end(x, 'stop', 'start'), axis=1)

    all_mutations.to_csv(output_folder + '/all_mutations_loc.csv',
                         sep=',',
                         index=False)

    print(all_mutations.head())
    print(all_mutations.shape)
    print(all_mutations.columns)

    df_complex = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name',
                                        'cancer']).agg({
                                                        'pos': ['nunique', 'count'],
                                                        'no_of_loc': concat_ints
                                                        }).reset_index()
    df_complex.columns = ['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name', 'cancer',
                          'pos_nunique', 'pos_count', 'no_of_loc']
    df_complex['complex'] = df_complex['pos_count'].apply(lambda x: 0 if x < 2 else 1)

    df_complex.to_csv(output_folder + '/complex.csv',
                      sep=',',
                      index=False)

    df_by_gene = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'cancer']).agg(
        {'indiv_name': ['nunique',
                        'count'],
         'pos': 'nunique',
         'no_of_loc': concat_ints
         }).reset_index()
    df_by_gene.columns = ['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'cancer', 'indiv_name_nunique',
                          'indiv_name_count', 'pos_nunique', 'no_of_loc']
    df_by_gene['if_complex'] = df_by_gene.apply(lambda x: if_complex(x, df_complex), axis=1)
    df_by_gene.to_csv(output_folder + '/occur.csv',
                      sep=',',
                      index=False)

    df_by_mutation = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'pos', 'ref', 'alt',
                                            'mutation_type', 'type_of_subst', 'cancer']).agg(
        {'indiv_name': 'nunique',
         'no_of_loc': concat_ints
         }).reset_index()

    df_by_mutation.columns = [' '.join(col).strip() for col in df_by_mutation.columns.values]

    df_by_mutation.to_csv(output_folder + '/by_distinct_mutation.csv',
                          sep=',',
                          index=False)
