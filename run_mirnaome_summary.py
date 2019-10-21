import pandas as pd
import os
import scipy.stats as ss
import click


@click.command()
@click.argument('output_folder')
def main(output_folder):
    files = [[(x[0] + '/' + y, x[0].split('/')[-1].replace('DATA_', '').replace('_22', '')) for y in x[2]]
             for x in os.walk(output_folder)
             if '_22' in x[0].split('/')[-1] or 'pancancer' in x[0].split('/')[-1]]
    flat_files = [file for sublist in files for file in sublist]
    csv_files = [file for file in flat_files if file[0].endswith('csv')]
    df = pd.DataFrame()

    for all_mutations_filtered_mut_type_gene, folder_name in [
        file for file in csv_files if file[0].endswith('all_mutations.csv')
    ]:
        df_temp = pd.read_csv(all_mutations_filtered_mut_type_gene)
        mutations_no = df_temp.groupby(['chrom', 'pos', 'indiv_name']).count().shape[0]
        localizations_no = df_temp.groupby(['chrom', 'pos']).count().shape[0]
        row = pd.DataFrame([[folder_name, df_temp.shape[0], mutations_no, localizations_no]],
                           columns=['folder', 'mutations_all_file_count', 'mutations_count',
                                    'mutation_positions'])
        df = pd.concat([df, row])
    df.set_index('folder', inplace=True)
    df.sort_index(inplace=True)

    df_complex = pd.DataFrame()
    for complex_file, folder_name in [
        file for file in csv_files if file[0].endswith('complex.csv')
    ]:
        df_temp = pd.read_csv(complex_file)
        df_temp = df_temp[df_temp['seq_type'] == 'mirna']
        df_temp.drop_duplicates(subset=['chrom', 'pre_name', 'indiv_name'], keep='first', inplace=True)

        row = pd.DataFrame([[folder_name, df_temp.shape[0]]],
                           columns=['folder', 'mutated_gene_patient_number'])
        df_complex = pd.concat([df_complex, row])
    df_complex.set_index('folder', inplace=True)
    df = df.join(df_complex)

    df_occur = pd.DataFrame()
    coords = pd.read_csv(
        '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates.bed',
        sep='\t', header=None)

    coords['length'] = coords[2] - coords[1]

    all_mirnas_len = coords[coords[3].str.contains('hsa-')]['length'].sum()

    coords['prob'] = coords['length'] / all_mirnas_len

    coords.to_excel(
        '/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates_length.xlsx',
        header=None)

    for occur, folder_name in [
        file for file in csv_files if file[0].endswith('occur.csv')
    ]:
        df_temp = pd.read_csv(occur)
        df_temp = df_temp[df_temp['seq_type'] == 'mirna']
        mutations = df.loc[folder_name, 'mutations_count']
        mutations_without_complex = df.loc[folder_name, 'mutated_gene_patient_number']

        if df_temp.shape[0] > 0:
            df_temp['binom'] = df_temp.apply(lambda x: round(ss.binom_test(x['indiv_name_count'],
                                                                     mutations,
                                                                     coords[coords[3] == x['pre_name']]['prob'].values[0],
                                                                     alternative='greater'), 8), axis=1)
            binom_values = list(df_temp['binom'].unique())
            binom_values.sort()
            df_temp['order'] = df_temp['binom'].apply(lambda x: binom_values.index(x) + 1)
            df_temp.sort_values('order', inplace=True)
            df_temp['FDR'] = df_temp.apply(lambda x: x['order'] / 1918 * 0.05, axis=1)
            df_temp['signif'] = df_temp['binom'] < df_temp['FDR']
            max_order = df_temp.loc[df_temp['signif'], 'order'].max()
            df_temp['BH'] = df_temp.apply(lambda x: x['binom'] * 1918 / x['order'], axis=1)
            mirnas = df_temp.loc[df_temp['order'] <= max_order, 'pre_name'].unique()

            df_temp['binom_complex'] = df_temp.apply(lambda x: round(ss.binom_test(x['indiv_name_nunique'],
                                                                             mutations_without_complex,
                                                                             coords[coords[3] == x['pre_name']]['prob'].values[0],
                                                                             alternative='greater'), 8), axis=1)
            binom_values = list(df_temp['binom_complex'].unique())
            binom_values.sort()
            df_temp['order_complex'] = df_temp['binom_complex'].apply(lambda x: binom_values.index(x) + 1)
            df_temp.sort_values('order_complex', inplace=True)
            df_temp['FDR_complex'] = df_temp.apply(lambda x: x['order_complex'] / 1918 * 0.05, axis=1)
            df_temp['signif_complex'] = df_temp['binom_complex'] < df_temp['FDR_complex']
            max_order = df_temp.loc[df_temp['signif_complex'], 'order_complex'].max()
            df_temp['BH_complex'] = df_temp.apply(lambda x: x['binom_complex'] * 1918 / x['order_complex'], axis=1)
            mirnas_complex = df_temp.loc[df_temp['order_complex'] <= max_order, 'pre_name'].unique()

            if folder_name == 'LUAD':
                df_temp.to_csv('test_luad.csv')

            if folder_name == 'LUSC':
                df_temp.to_csv('test_lusc.csv')

        else:
            mirnas = []
            mirnas_complex = []
        mirnas_count = len(mirnas)
        mirnas_complex_count = len(mirnas_complex)
        row = pd.DataFrame([[folder_name, df_temp.shape[0], mirnas, mirnas_count,
                             mirnas_complex, mirnas_complex_count]], columns=['folder', 'mutated_genes_count',
                                                                              'significant_mirnas',
                                                                              'significant_mirnas_count',
                                                                              'significant_mirnas_complex',
                                                                              'significant_mirnas_complex_count'])
        df_occur = pd.concat([df_occur, row])
    df_occur.set_index('folder', inplace=True)
    df = df.join(df_occur)
    if output_folder.endswith('/'):
        df.to_excel(output_folder + 'summary.xlsx')
    else:
        df.to_excel(output_folder + '/summary.xlsx')


if __name__ == "__main__":
    main()
