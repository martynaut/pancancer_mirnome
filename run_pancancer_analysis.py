import pandas as pd
import click
import os


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
def main(input_folder, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    files = [[(x[0] + '/' + y, x[0].split('/')[-1].replace('DATA_', '').replace('_22', '')) for y in x[2]] for x in
             os.walk(input_folder)
             if '_22' in x[0].split('/')[-1]]
    flat_files = [file for sublist in files for file in sublist]
    csv_files = [file for file in flat_files if file[0].endswith('csv') and (
                 'pancancer' not in file[1])]

    df_fiels_summary = pd.DataFrame()
    for summary_file, folder_name in [
        file for file in csv_files if file[0].endswith('files_summary_count_per_patient.csv')
    ]:
        df_temp = pd.read_csv(summary_file)
        df_temp = df_temp.groupby('indiv_name').count()
        row = pd.DataFrame([[folder_name, df_temp.shape[0]]], columns=['folder', 'number_of_patients'])
        df_fiels_summary = pd.concat([df_fiels_summary, row], sort=False)
    if output_folder.endswith('/'):
        df_fiels_summary.to_excel(output_folder + 'number_of_patients.xlsx', index=False)
    else:
        df_fiels_summary.to_excel(output_folder + '/number_of_patients.xlsx', index=False)

    df_count_summary = pd.DataFrame()
    for count_file, folder_name in [
        file for file in csv_files if file[0].endswith('patient_mutation_count.csv')
    ]:
        df_temp = pd.read_csv(count_file)
        df_temp['cancer_id'] = folder_name
        df_count_summary = pd.concat([df_count_summary, df_temp], sort=False)
    if output_folder.endswith('/'):
        df_count_summary.to_excel(output_folder + 'patient_mutation_count.xlsx', index=False)
    else:
        df_count_summary.to_excel(output_folder + '/patient_mutation_count.xlsx', index=False)

    all_mutation_count = pd.DataFrame()
    for count_file, folder_name in [
        file for file in csv_files if file[0].endswith('patient_mutation_count.csv')
    ]:
        df_temp = pd.read_csv(count_file)
        df_temp['cancer'] = folder_name
        all_mutation_count = pd.concat([all_mutation_count, df_temp], sort=False)
    if output_folder.endswith('/'):
        all_mutation_count.to_excel(output_folder + 'number_of_mutations_all_patients.xlsx', index=False)
    else:
        all_mutation_count.to_excel(output_folder + '/number_of_mutations_all_patients.xlsx', index=False)

    df_all_mutation = pd.DataFrame()
    for all_mutations_filtered_mut_type_gene, folder_name in [
        file for file in csv_files if file[0].endswith('all_mutations.csv')
    ]:
        print(folder_name)
        df_temp = pd.read_csv(all_mutations_filtered_mut_type_gene)
        df_temp['cancer'] = folder_name
        df_all_mutation = pd.concat([df_all_mutation, df_temp], sort=False)
    df_all_mutation.to_csv(output_folder + 'all_mutations.csv')

    df_complex = pd.DataFrame()
    for complex_file, folder_name in [
        file for file in csv_files if file[0].endswith('complex.csv')
    ]:
        df_temp = pd.read_csv(complex_file)
        df_temp['cancer'] = folder_name
        df_complex = pd.concat([df_complex, df_temp], sort=False)
    df_complex.to_csv(output_folder + 'complex.csv')

    df_occur = pd.DataFrame()
    for occur_file, folder_name in [
        file for file in csv_files if file[0].endswith('occur.csv')
    ]:
        df_temp = pd.read_csv(occur_file)
        df_occur = pd.concat([df_occur, df_temp], sort=False)
    df_occur = df_occur.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type'], as_index=False).agg({
        'indiv_name_nunique': sum,
        'indiv_name_count': sum,
        'pos_nunique': sum,
        'if_complex': lambda x: 1 if sum(x) > 0 else 0
    })
    df_occur.to_csv(output_folder + 'occur.csv')


if __name__ == "__main__":
    main()
