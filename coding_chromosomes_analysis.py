import gzip
import pandas as pd
import os
from prepare_vcf_files_helpers import update_dict_with_file


pd.options.mode.chained_assignment = None


def each_file_processing_coding(filename, reference):

    print(filename)

    dataframe_records = pd.DataFrame()

    with gzip.open(filename) as f:
        columns = []
        for line in f.readlines():
            line = line.decode('ascii')
            if line[:1] == '##':
                pass
            elif line[:1] == '#':
                columns = line.replace('#', '').strip().split('\t')
            else:
                position = line.split('\t')[:5]
                if reference[(reference['chr'] == position[0]) &
                             (reference['start_ref'] < int(position[1])) &
                             (reference['stop_ref'] > int(position[1]))].shape[0] > 0:

                    new_record = pd.DataFrame([line.replace('\n',
                                                            '').replace(';',
                                                                        ':').replace('"',
                                                                                     '').split('\t')],
                                              columns=columns)
                    if new_record['FILTER'].str.contains('PASS').any():
                        new_record['chr'] = new_record['CHROM']
                        new_record['count'] = 1
                        new_record = new_record[['chr', 'count']]
                        dataframe_records = pd.concat([dataframe_records, new_record])
    try:
        dataframe_records = dataframe_records.groupby('chr', as_index=False).sum()
    except KeyError:
        new_record = pd.DataFrame([['chrX', 0]], columns=['chr', 'count'])
        dataframe_records = pd.concat([dataframe_records, new_record])
        dataframe_records = dataframe_records.groupby('chr', as_index=False).sum()
    return dataframe_records


def all_files_processing_coding(input_folder, output_folder, coordinates_file):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    coordinates = pd.read_table(coordinates_file,
                                names=['chr', 'start', 'stop', 'gene'])

    coordinates['start_ref'] = coordinates['start']
    coordinates['stop_ref'] = coordinates['stop']

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
    files = [file for sublist in files for file in sublist]
    files_temp = [[x[0] + '/' + y for y in x[2]] for x in os.walk(output_folder + '/temp')]
    files_temp = [file for sublist in files_temp for file in sublist]

    files.extend(files_temp)

    with open(output_folder + '/do_not_use.txt', 'r') as file_dont:
        do_not_use = []
        for line in file_dont.readlines():
            do_not_use.append(line[:-1])
        gz_files = []
        for file in files:
            if file[-6:] == 'vcf.gz' and file not in do_not_use:
                gz_files.append(file)
            else:
                pass

    # files summary
    dict_with_files = {}
    for gz_file in gz_files:
        dict_with_files = update_dict_with_file(gz_file, dict_with_files)

    df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')

    df_dict_with_files.index.name = 'filename'
    df_dict_with_files.to_csv(output_folder + '/files_summary.csv', sep=',')

    df_dict_with_files_grouped = df_dict_with_files.reset_index().groupby(['indiv_name',
                                                                           'type_of_file']).agg('nunique')
    df_dict_with_files_grouped.to_csv(output_folder + '/files_summary_count_per_patient.csv', sep=',')

    df_file_type = df_dict_with_files.reset_index().groupby(['type_of_file'])[['filename']].agg('nunique')
    df_file_type.columns = ['count_filename']
    df_file_type.to_csv(output_folder + '/files_count_per_type.csv', sep=',')

    counter = 0
    all_files = df_dict_with_files.shape[0]
    for file_type in list(df_dict_with_files['type_of_file'].unique()):

        results_df = pd.DataFrame()
        for x in range(1, 23):
            new_row = pd.DataFrame([['chr' + str(x), 0]], columns=['chr', 'zeros'])
            results_df = pd.concat([results_df, new_row])
        new_row = pd.DataFrame([['chrX', 0]], columns=['chr', 'zeros'])
        results_df = pd.concat([results_df, new_row])
        new_row = pd.DataFrame([['chrY', 0]], columns=['chr', 'zeros'])
        results_df = pd.concat([results_df, new_row])
        for file in list(df_dict_with_files.loc[df_dict_with_files['type_of_file'] == file_type, :].index):
            counter += 1
            print(str(counter) + ' / ' + str(all_files) + '\n')
            dataframe = each_file_processing_coding(file,
                                                    coordinates)
            results_df = results_df.join(dataframe.set_index('chr'), on=['chr'], rsuffix=file)

        results_df.to_csv(output_folder + '/results_{}.csv'.format(file_type),
                          sep=',',
                          index=False)


def merge_data(output_folder, coordinates_file):
    coordinates = pd.read_table(coordinates_file,
                                names=['chr', 'start', 'stop', 'gene'])
    coordinates.drop_duplicates(subset=['chr', 'start', 'stop'], keep='first', inplace=True)
    coordinates['length'] = coordinates['stop'] - coordinates['start']
    df_joined = coordinates.groupby('chr', as_index=False).sum()[['chr', 'length']]

    df_varscan = pd.read_csv(output_folder + '/results_varscan2.csv')
    df_muse = pd.read_csv(output_folder + '/results_muse.csv')
    df_ss = pd.read_csv(output_folder + '/results_somaticsniper.csv')
    df_mutect = pd.read_csv(output_folder + '/results_mutect2.csv')

    df_varscan = pd.DataFrame(df_varscan.fillna(0).set_index('chr').sum(axis=1), columns=['varscan'])
    df_muse = pd.DataFrame(df_muse.fillna(0).set_index('chr').sum(axis=1), columns=['muse'])
    df_ss = pd.DataFrame(df_ss.fillna(0).set_index('chr').sum(axis=1), columns=['ss'])
    df_mutect = pd.DataFrame(df_mutect.fillna(0).set_index('chr').sum(axis=1), columns=['mutect'])
    df_joined = df_joined.join(df_varscan, on='chr').join(df_muse, on='chr').join(df_ss, on='chr')
    df_joined = df_joined.join(df_mutect, on='chr')
    for type in ['varscan', 'muse',
                 'ss',
                 'mutect'
                 ]:
        df_joined[type + '_freq'] = df_joined[type] / df_joined['length'] * 1000
    print(df_joined)
    df_joined.to_csv(output_folder + '/chromosomes_summary.csv')


