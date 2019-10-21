import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/all_mutations_weights.csv'
)

print(df.head())

no_of_cancers = df['cancer'].nunique()

df_accepted = df[df['>10000'] != 1]
df_rejected = df[df['>10000'] == 1]

df_grouped_mutations = df.groupby(['cancer', 'chrom',
                                   'pos', 'indiv_name'], as_index=False).count()

df_grouped_cancer = df_grouped_mutations.groupby('cancer', as_index=False).agg(
    {'indiv_name': 'nunique',
     'pos': 'count'
     }
)

df_grouped_cancer['freq'] = df_grouped_cancer['pos'] / df_grouped_cancer['indiv_name']

df_grouped_patient = df_grouped_mutations.groupby(['cancer', 'indiv_name'], as_index=False).count()

fig, axs = plt.subplots(3, no_of_cancers, sharey='row', figsize=(30, 20))


n = 0
for cancer in df['cancer'].unique():
    df_temp = df_grouped_patient[df_grouped_patient['cancer'] == cancer].copy()
    df_temp.sort_values('pos', inplace=True)
    df_temp.reset_index(inplace=True)
    axs[0][n].plot(df_temp.index, df_temp['pos'], linestyle='None', marker='.')
    axs[0][n].set_xlabel(cancer + '\n(n={})'.format(df_temp['indiv_name'].nunique()))
    n = n + 1

axs[0][0].set_ylabel('Number of mutations in miRNA genes per Exome')
axs[1][0].set_ylabel('Total number of mutations per Exome')
axs[2][0].set_ylabel('Ratio of miRNA genes mutations\nand total number of mutations per Exome')

# plt.yscale('log')

plt.savefig(
    '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/miRNA_pancancer_mutation_count_10000.png',
    type='png', bbox_inches='tight')
plt.show()

fig, axs = plt.subplots(1, 1, figsize=(4, 10))

axs.barh(range(no_of_cancers), df_grouped_cancer['freq'],
         tick_label=df_grouped_cancer['cancer'])
axs.set_xlabel('Average number of mutations\nin miRNome per Exome')

plt.show()
