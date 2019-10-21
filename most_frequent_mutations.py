import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

SMALL_SIZE = 7
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title



df = pd.read_csv('/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/all_mutations.csv')

# df_grouped = df.groupby(['chrom', 'pos', 'indiv_name', 'pre_name'], as_index=False).count()
#
# # mutations
# df_count = df_grouped[['chrom', 'pos', 'indiv_name', 'pre_name']].groupby(['chrom', 'pos',
#                                                                            'pre_name'],
#                                                                           as_index=False).count()
# df_count.sort_values('indiv_name', ascending=False, inplace=True)
# df_head = df_count.head(20).copy()
# df_head['mut_name'] = df_head.apply(lambda x: x['pre_name'] + ' ' + x['chrom'] + ':' + str(x['pos']),
#                                     axis=1)
#
# # mirnas
#
# df_count_mi = df_grouped[['chrom', 'pos', 'indiv_name', 'pre_name']].groupby(['chrom', 'pre_name'],
#                                                                              as_index=False).agg({
#     'indiv_name': 'nunique'
# })
# df_count_mi.sort_values('indiv_name', ascending=False, inplace=True)
# df_head_mi = df_count_mi.head(20).copy()
#
#
# fig, axs = plt.subplots(1, 2, figsize=(14, 10))
#
# axs[0].barh(range(1, 21), df_head['indiv_name'],
#          tick_label=df_head['mut_name'])
# axs[0].set_ylim(20.5, 0.5)
# axs[0].set_title('Number of occurances\nof most frequent mutations in all cancer types')
#
# axs[1].barh(range(1, 21), df_head_mi['indiv_name'],
#          tick_label=df_head_mi['pre_name'])
# axs[1].set_ylim(20.5, 0.5)
# axs[1].set_title('Number of patients with mutation\nin most frequently mutated miRNA in all cancer types')
#
# plt.tight_layout()
#
# plt.savefig(
#         '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/most_freq_mutations.png',
#         type='png', bbox_inches='tight')
# plt.savefig(
#         '/Users/martynaurbanek/Documents/private/repos/files/output/DATA_pancancer/most_freq_mutations.tiff',
#         format='tiff', dpi=300, transparent=True, pil_kwargs={'compression': "jpeg"})


df_grouped2 = df.groupby(['chrom', 'pos', 'indiv_name', 'pre_name', 'cancer'], as_index=False).count()
df_count2 = df_grouped2[['chrom', 'pos', 'indiv_name', 'pre_name', 'cancer']].groupby(['chrom',
                                                                           'pre_name', 'cancer'],
                                                                          as_index=False).agg({
    'indiv_name': 'nunique'
})

df_index = df_count2.groupby(['chrom', 'pre_name'], as_index=False).sum()
df_index.sort_values('indiv_name', ascending=False, inplace=True)
df_index = df_index[['chrom', 'pre_name']].head(20)

df_count2 = df_index.join(df_count2.set_index(['chrom', 'pre_name']),
                           ['chrom', 'pre_name'],
                           'inner')


print(df_count2.head())
plot_dw = df_count2[['pre_name', 'cancer', 'indiv_name']].set_index(['pre_name', 'cancer'])
plot_dw = plot_dw.unstack()
plot_dw = plot_dw['indiv_name']
print(plot_dw)
plot_dw_index = plot_dw.sum(axis=1)
plot_dw['sorter'] = plot_dw_index
plot_dw.sort_values('sorter', ascending=True, inplace=True)
plot_dw.drop('sorter', axis=1, inplace=True)

colors = plt.cm.tab20b.colors
print(colors)
colors2 = plt.cm.tab20c.colors
colors = colors + colors2
cmap = matplotlib.colors.ListedColormap(name='mycmap', colors=colors)


plot_dw.plot(kind='barh', stacked=True, colormap=cmap)
plt.legend(bbox_to_anchor=(0.9, -0.1), ncol=6)
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()


df_grouped2 = df.groupby(['chrom', 'pos', 'indiv_name', 'pre_name', 'cancer'], as_index=False).count()
df_grouped2['mut_name'] = df_grouped2.apply(lambda x: x['pre_name'] + ' ' + x['chrom'] + ':' + str(x['pos']),
                                    axis=1)
df_count2 = df_grouped2[['mut_name', 'indiv_name', 'cancer']].groupby(['mut_name', 'cancer'],
                                                                          as_index=False).agg({
    'indiv_name': 'nunique'
})

df_index = df_count2.groupby(['mut_name'], as_index=False).sum()
df_index.sort_values('indiv_name', ascending=False, inplace=True)
df_index = df_index[['mut_name']].head(20)

df_count2 = df_index.join(df_count2.set_index(['mut_name']),
                           ['mut_name'],
                           'inner')


print(df_count2.head())
plot_dw = df_count2[['mut_name', 'cancer', 'indiv_name']].set_index(['mut_name', 'cancer'])
plot_dw = plot_dw.unstack()
plot_dw = plot_dw['indiv_name']
print(plot_dw)
plot_dw_index = plot_dw.sum(axis=1)
plot_dw['sorter'] = plot_dw_index
plot_dw.sort_values('sorter', ascending=True, inplace=True)
plot_dw.drop('sorter', axis=1, inplace=True)

colors = plt.cm.tab20b.colors
print(colors)
colors2 = plt.cm.tab20c.colors
colors = colors + colors2
cmap = matplotlib.colors.ListedColormap(name='mycmap', colors=colors)


plot_dw.plot(kind='barh', stacked=True, colormap=cmap)
plt.legend(bbox_to_anchor=(0.9, -0.1), ncol=6)
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()

