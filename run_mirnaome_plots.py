import pandas as pd
import matplotlib
import click
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import os


@click.command()
@click.option('--output_folder', '-o')
def main(output_folder='/Users/martynaurbanek/Documents/private/repos/files/output/'):
    if output_folder is None:
        output_folder = '/Users/martynaurbanek/Documents/private/repos/files/output/'

    files = [[(x[0] + '/' + y, x[0].split('/')[-1].replace('DATA_', '')) for y in x[2]] for x in os.walk(output_folder)]
    flat_files = [file for sublist in files for file in sublist]
    csv_files = [file for file in flat_files if file[0].endswith('csv')]

    folders = [
            'ACC_test',
            'BLCA_test',
            'BRCA_test',
            'CESC_test',
            'CHOL_test',
            'COAD_test',
            'COAD2_test',
            'DLBC_test',
            'ESCA_test',
            'GBM_test',
            'KICH_test',
            'KIRC_test',
            'KIRP_test',
            'LAML_test',
            'LGG_test',
            'LIHC_test',
             'LUAD_22',
            'LUSC_new',
            'MESO_test',
            'OV_test',
            'PAAD_test',
            'PCPG1_test',
            'PCPG2_test',
            'PRAD_test',
            'READ1_test2',
            'READ2_test',
            'SARC1_test',
            'SARC2_test',
            'SARC3_test',
            'SKCM_test',
            'STAD_test',
            'TCGT_test',
            'THCA_test',
            'THYM1_test',
            'THYM2_test',
            'UCEC_test',
            'UCS_test',
            'UVM_test',
            'pancancer'
        ]

    df = pd.read_excel(output_folder + 'summary.xlsx')
    df = df[df['folder'].isin(
        folders
    )]
    df = df[['folder', 'significant_mirnas']]
    df['mirnas'] = df['significant_mirnas'].apply(
        lambda x: [y.strip() for y in x.replace('[', '').replace(']', '').
                   replace("'", "").replace("'", "").split(' ')]
    )
    mirnas = list(set([mirna for sublist in df['mirnas'].values for mirna in sublist]))
    try:
        mirnas.remove('')
    except ValueError:
        pass
    for mirna in mirnas:
        df[mirna] = df['mirnas'].apply(lambda x: 1 if mirna in x else 0)
    df.set_index('folder', inplace=True)
    df = df[mirnas]
    plt.figure(figsize=(50, 10))
    sns.heatmap(df, xticklabels=True, yticklabels=True, cbar=False,
                cmap='Greys')
    matplotlib.pyplot.savefig(output_folder + 'plots/plot_heatmap.png', dpi=300, bbox_inches='tight')

    cm = sns.clustermap(df, xticklabels=True, yticklabels=True,
                        cmap='Greys',
                        row_cluster=True,
                        col_cluster=False,
                        figsize=(50, 10),
                        method='ward')
    cm.cax.set_visible(False)
    matplotlib.pyplot.savefig(output_folder + 'plots/plot_clustermap_ward.png', dpi=300, bbox_inches='tight',
                              method='ward')


if __name__ == "__main__":
    main()
