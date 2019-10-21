import pandas as pd
import matplotlib
import click
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns


@click.command()
@click.option('--output_folder', '-o')
def main(output_folder='/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/'):
    if output_folder is None:
        output_folder = '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/'

    df_count = pd.read_excel(output_folder + 'number_of_patients.xlsx')

    df = pd.read_csv(output_folder + 'all_mutations.csv')
    df = df.join(df_count.set_index('folder'), on='cancer', how='left')
    # print(df.head())
    df = df.groupby(['cancer', 'pre_name', 'number_of_patients'], as_index=False).agg({'indiv_name': 'nunique'})
    # print(df.head())

    # df_temp = df.groupby('pre_name', as_index=False).agg({'indiv_name': 'nunique'})
    # df_temp['number_of_patients'] = df_count['number_of_patients'].sum()
    # df_temp['cancer'] = 'pancancer'
    # print(df_temp.head())

    # df = pd.concat([df, df_temp], sort=True)

    df['freq'] = df['indiv_name'] / df['number_of_patients'] * 100

    df_pivot = df.pivot(index='cancer', columns='pre_name', values='freq')
    df_pivot.fillna(0, inplace=True)
    df_pivot.to_excel(output_folder + 'pivot_mutation_freq.xlsx')
    # print(df_pivot.tail())
    plt.figure(figsize=(50, 10))
    sns.heatmap(df_pivot, yticklabels=True, cmap="PuRd")
    matplotlib.pyplot.savefig(output_folder + 'plots/plot_heatmap_freq_all.png', dpi=300, bbox_inches='tight')

    df2 = pd.read_excel(output_folder + '../summary.xlsx')

    df2 = df2[['folder', 'significant_mirnas']]
    df2['mirnas'] = df2['significant_mirnas'].apply(
        lambda x: [y.strip() for y in x.replace('[', '').replace(']', '').
                   replace("'", "").replace("'", "").split(' ')]
    )
    mirnas = list(set([mirna for sublist in df2['mirnas'].values for mirna in sublist]))
    try:
        mirnas.remove('')
    except ValueError:
        pass

    df_pivot2 = df_pivot[mirnas]
    plt.figure(figsize=(50, 10))
    sns.heatmap(df_pivot2, xticklabels=True, yticklabels=True, cmap="Blues")
    matplotlib.pyplot.savefig(output_folder + 'plots/plot_heatmap_freq_significant.png', dpi=300, bbox_inches='tight')

    pre_mirna_count = df[['pre_name', 'indiv_name']].groupby('pre_name',
                                                             as_index=False).sum()

    mirnas2 = pre_mirna_count[pre_mirna_count['indiv_name'] > 15]['pre_name'].unique()
    print(mirnas2)
    print(len(mirnas2))


    df_pivot3 = df_pivot[mirnas2]
    plt.figure(figsize=(50, 10))
    sns.heatmap(df_pivot3, xticklabels=True, yticklabels=True, cmap="Blues")
    plt.ylabel('')
    plt.xlabel('')
    matplotlib.pyplot.savefig(output_folder + 'plots/plot_heatmap_frequent.png', dpi=300, bbox_inches='tight',
                              type='png')
    matplotlib.pyplot.savefig(output_folder + 'plots/plot_heatmap_frequent.jpeg', dpi=300, bbox_inches='tight',
                              type='jpeg')


if __name__ == "__main__":
    main()
