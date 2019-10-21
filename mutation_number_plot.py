import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def main():

    df_all_mut = pd.read_excel(
        '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/' +
        'number_of_mutations_all_patients.xlsx')
    print(df_all_mut.head())

    df = pd.read_csv('/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/all_mutations.csv')

    print(df.head())

    no_of_cancers = df['cancer'].nunique()

    df_grouped_mutations = df.groupby(['cancer', 'chrom',
                                       'pos', 'indiv_name'],
                                      as_index=False).count()

    df_all_mut_joined = df_all_mut.join(df_grouped_mutations[['indiv_name',
                                                              'cancer', 'pos']].groupby(['indiv_name',
                                                                                         'cancer']).count(),
                                        ['patient_id', 'cancer'], 'outer')
    df_all_mut_joined['pos'] = df_all_mut_joined['pos'].fillna(0)
    df_all_mut_joined['ratio'] = df_all_mut_joined['pos'] / df_all_mut_joined['mutation_count']
    df_all_mut_joined['color'] = df_all_mut_joined['mutation_count'].apply(
        lambda x: 'red' if x > 10000 else 'blue'
    )
    df_all_mut_joined.to_csv(
        '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/patient_mutations.csv',
        index=False)

    df_grouped_cancer = df_grouped_mutations.groupby('cancer', as_index=False).agg(
        {'indiv_name': 'nunique',
         'pos': 'count'
         }
    )

    df_grouped_cancer['freq'] = df_grouped_cancer['pos'] / df_grouped_cancer['indiv_name']

    # fig, axs = plt.subplots(3, no_of_cancers, sharey='row', figsize=(30, 20))

    # n = 0
    # for cancer in df.sort_values('cancer')['cancer'].unique():
    #     df_temp = df_all_mut_joined[df_all_mut_joined['cancer'] == cancer].copy()
    #     df_temp.sort_values(['pos'], inplace=True)
    #     df_temp.reset_index(inplace=True)
    #
    #     axs[0][n].scatter(df_temp.index, df_temp['pos'], color=df_temp['color'])
    #     axs[0][n].set_xlabel(cancer + '\n(n={})'.format(df_temp.shape[0]))
    #     axs[0][n].set_yscale('symlog', linthreshy=10)
    #
    #     df_temp.sort_values(['mutation_count'], inplace=True)
    #     df_temp.drop('index', axis=1, inplace=True)
    #     df_temp.reset_index(inplace=True)
    #
    #     axs[1][n].scatter(df_temp.index, df_temp['mutation_count'], color=df_temp['color'])
    #     axs[1][n].set_xlabel(cancer + '\n(n={})'.format(df_temp.shape[0]))
    #     axs[1][n].set_yscale('symlog', linthreshy=10)
    #
    #     df_temp.sort_values(['ratio'], inplace=True)
    #     df_temp.drop('index', axis=1, inplace=True)
    #     df_temp.reset_index(inplace=True)
    #
    #     axs[2][n].scatter(df_temp.index, df_temp['ratio'], color=df_temp['color'])
    #     axs[2][n].set_xlabel(cancer + '\n(n={})'.format(df_temp.shape[0]))
    #     axs[2][n].set_yscale('symlog', linthreshy=0.01)
    #     n = n + 1
    # fig.align_ylabels(axs)
    # axs[0][0].set_ylabel('Number of mutations in miRNA genes per Exome')
    # axs[1][0].set_ylabel('Total number of mutations per Exome')
    # axs[2][0].set_ylabel('Ratio of miRNA genes mutations\nand total number of mutations per Exome')
    #
    # # plt.yscale('log')
    #
    # plt.savefig(
    #     '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/miRNA_pancancer_mutation_count.png',
    #     type='png', bbox_inches='tight')
    # plt.savefig(
    #     '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/miRNA_pancancer_mutation_count.tiff',
    #     format='tiff', dpi=300, transparent=True, pil_kwargs={'compression': "jpeg"})
    # plt.show()

    fig, axs = plt.subplots(6, 7, figsize=(12, 12))

    gs = axs[0, 0].get_gridspec()
    # remove the underlying axes
    for ax in axs[:2, :2]:
        for axx in ax:
            axx.remove()
    axbig = fig.add_subplot(gs[:2, :2])

    axbig.scatter(df_all_mut_joined['mutation_count'],
                  df_all_mut_joined['pos'],
                  color=df_all_mut_joined['color'])
    axbig.set_xscale('symlog', linthreshx=10, linscalex=0.5)
    # axbig.set_yscale('symlog', linthreshy=10, linscaley=0.1)

    axbig.set_xlabel('pancancer' + '\n(n={})'.format(df_all_mut_joined.shape[0]))

    n = 2
    m = 0

    for cancer in df.sort_values('cancer')['cancer'].unique():
        df_temp = df_all_mut_joined[df_all_mut_joined['cancer'] == cancer].copy()
        axs[n][m].scatter(df_temp['mutation_count'], df_temp['pos'],
                          color=df_temp['color'])
        axs[n][m].set_xlabel(cancer + '\n(n={})'.format(df_temp.shape[0]))
        axs[n][m].set_xscale('symlog', linthreshx=10, linscalex=0.5)
        axs[n][m].set_yticks(np.arange(0, df_temp['pos'].max()+1,
                             step=round(df_temp['pos'].max()/4)+1))
        n = n+1
        if n == 6:
            if m > 0:
                n = 0
                m = m+1
            else:
                n = 2
                m = m+1


    plt.show()


if __name__ == "__main__":
    main()
