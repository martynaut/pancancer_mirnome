import pandas as pd
import matplotlib
import click
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from matplotlib import lines


def create_plot(data_df, output_name, mutations=0, genes=0):

    df_plot_5 = prepare_data_5p(data_df)

    df_plot_3 = prepare_data_3p(data_df)

    max_value = max(df_plot_5['pos'].max(), df_plot_3['pos'].max())

    rc("pdf", fonttype=42)
    sns.set_style(style='white')

    palette = {'flanking-5': 'grey',
               'flanking-3': 'grey',
               'pre-seed': 'lightblue',
               'seed': 'darkblue',
               'post-seed': 'steelblue',
               'loop': 'grey'
               }

    with PdfPages(output_name) as pdf:

        fig = plt.figure(figsize=(25, 10))
        ax = fig.add_axes([0, 0, 1, 1])

        ax.axis('off')
        im = plt.imread('/Users/martynaurbanek/Documents/private/repos/files/primirna_background.tiff')
        plt.imshow(im)
        plt.xticks([])
        plt.yticks([])

        y_min, y_max = ax.get_ylim()
        x_min, x_max = ax.get_xlim()

        plt.text(x_max * 0.98, y_min * 1.07, str(mutations), horizontalalignment='center',
                 verticalalignment='center', fontdict={'size': '19'})
        plt.text(x_max * 0.98, y_min * 1.14, str(genes), horizontalalignment='center',
                 verticalalignment='center', fontdict={'size': '19'})
        line = lines.Line2D([x_max * 0.97, x_max * 0.99], [y_min * 1.105, y_min * 1.105], lw=1.5, color='black')
        ax.add_line(line)
        line.set_clip_on(False)

        plt.text(x_max * 0.26, y_min * -0.09, 'flanking region', horizontalalignment='center', fontdict={'size': '19'})
        line = lines.Line2D([x_max * 0.02, x_max * 0.5], [y_min * -0.07, y_min * -0.07], lw=1.5, color='grey', alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)

        plt.text(x_max * 0.70, y_min * -0.09, 'miRNA', horizontalalignment='center', fontdict={'size': '19'})
        line = lines.Line2D([x_max * 0.51, x_max * 0.89], [y_min * -0.07, y_min * -0.07], lw=1.5, color='steelblue', alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)

        plt.text(x_max * 0.94, y_min * -0.09, 'loop', horizontalalignment='center', fontdict={'size': '19'})
        line = lines.Line2D([x_max * 0.9, x_max * 0.98], [y_min * -0.07, y_min * -0.07], lw=1.5, color='grey', alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)

        plt.text(x_max * 0.585, y_min * -0.03, 'seed', horizontalalignment='center', fontdict={'size': '19'})
        line = lines.Line2D([x_max * 0.53, x_max * 0.64], [y_min * -0.01, y_min * -0.01], lw=1.5, color='darkblue', alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)

        line = lines.Line2D([x_max * 0.735, x_max * 0.835], [y_min * 1.125, y_min * 1.125], lw=1.5, color='darkblue',
                            alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)
        plt.text(x_max * 0.785, y_min * 1.175, 'seed', horizontalalignment='center', fontdict={'size': '19'})

        a = plt.axes([.058, .52, .93, .35])

        hue_order = ['flanking-5', 'pre-seed', 'seed', 'post-seed', 'loop']
        labels = [str(x) for x in range(-25, 0)] + \
                  [str(x) for x in range(1, 23)] + ['+' + str(x) for x in range(1, 6)]

        ax = sns.barplot(x="from_start", y="pos", hue="type",
                         data=df_plot_5, dodge=False,
                         hue_order=hue_order,
                         palette=palette,
                         ax=a)

        for loc in ['right', 'top', 'left', 'bottom']:
            ax.spines[loc].set_visible(False)

        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.legend(loc='upper right', ncol=2)

        plt.setp(ax.get_legend().get_texts(), fontsize='18')
        plt.setp(ax.get_legend().get_title(), fontsize='22')
        plt.grid(b=True, which='major', axis='y', color='lightgrey',
                 linestyle='-', linewidth=0.75, zorder=2, alpha=0.5)

        ax.set_xticklabels(labels)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(19)

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(19)


        if max_value > 5:
            plt.yticks(np.arange(0, max_value + 1, np.floor(max_value/3)))
            plot_limit = max_value + 2
        else:
            plt.yticks(np.arange(0, max_value + 1, 1))
            plot_limit = max_value + 1

        ax.set_ylim([0, plot_limit])
        ax.xaxis.tick_bottom()
        for label in ax.get_xticklabels():
            if label.get_text() not in ['-25', '-20', '-15', '-10', '-5', '1', '5', '10', '15', '20',
                                        '+1', '+5']:
                label.set_visible(False)

        new_width = 0.35
        for patch in ax.patches:
            current_width = patch.get_width()
            diff = current_width - new_width

            # we change the bar width
            patch.set_width(new_width)

            # we recenter the bar
            patch.set_x(patch.get_x() + diff * .5)

        ax.tick_params(axis='both', which='major', pad=8, width=0.5)

        plt.setp(ax.patches, linewidth=0)

        ax.get_legend().remove()

        b = plt.axes([.022, .02, .93, .35], facecolor='w')
        hue_order = ['flanking-3', 'pre-seed', 'seed', 'post-seed', 'loop']
        labels = ['+' + str(x) for x in range(1, 26)][::-1] + [str(x) for x in range(1, 23)][::-1] + \
                     ['-' + str(x) for x in range(1, 6)]

        ax = sns.barplot(x="from_start", y="pos", hue="type",
                         data=df_plot_3, dodge=False,
                         hue_order=hue_order,
                         palette=palette,
                         ax=b)

        for loc in ['right', 'top', 'left', 'bottom']:
            ax.spines[loc].set_visible(False)

        ax.set_xlabel('')
        ax.set_ylabel('')

        plt.grid(b=True, which='major', axis='y', color='lightgrey',
                 linestyle='-', linewidth=0.75, zorder=1, alpha=0.5)

        ax.set_xticklabels(labels)

        if max_value > 5:
            plt.yticks(np.arange(0, max_value + 1, np.floor(max_value/3)))
            plot_limit = max_value + 2
        else:
            plt.yticks(np.arange(0, max_value + 1, 1))
            plot_limit = max_value + 1

        ax.get_legend().remove()
        ax.set_ylim([plot_limit, 0])

        ax.xaxis.tick_top()

        for tick in ax.xaxis.get_major_ticks():
            tick.label2.set_fontsize(19)

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(19)

        for label in ax.get_xticklabels():
            if label.get_text() not in ['+25', '+20', '+15', '+10', '+5', '+1', '1', '5', '10', '15', '20',
                                        '-5']:
                label.set_visible(False)

        new_width = 0.35
        for patch in ax.patches:
            current_width = patch.get_width()
            diff = current_width - new_width

            # we change the bar width
            patch.set_width(new_width)

            # we recenter the bar
            patch.set_x(patch.get_x() + diff * .5)

        plt.setp(ax.patches, linewidth=0)

        plt.savefig(pdf, format='pdf', transparent=True)
        plt.close()


def prepare_data_5p(df_temp):
    dataframe = df_temp[(df_temp['arm'] == '5p') |
                        ((df_temp['arm'] == 'loop') & (df_temp['from_start'] < 6))].groupby(['arm', 'type',
                                                                                             'from_start'],
                                                                                            as_index=False
                                                                                            )[['pos']]. \
        count()[[
            'arm', 'type',
            'from_start',
            'pos'
        ]]

    dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '5p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '5p'), 'from_start'] + 25

    dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '5p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '5p'), 'from_start'] + 26

    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '5p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '5p'), 'from_start'] + 33

    dataframe.loc[(dataframe['type'] == 'loop') & (dataframe['arm'] == 'loop'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'loop') & (dataframe['arm'] == 'loop'), 'from_start'] + 47

    for x in range(1, 26):
        if dataframe[(dataframe['type'] == 'flanking-5') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'flanking-5', x, 0]], columns=['arm', 'type',
                                                                          'from_start',
                                                                          'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(26, 27):
        if dataframe[(dataframe['type'] == 'pre-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'pre-seed', x, 0]], columns=['arm', 'type',
                                                                        'from_start',
                                                                        'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(27, 34):
        if dataframe[(dataframe['type'] == 'seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'seed', x, 0]], columns=['arm', 'type',
                                                                    'from_start',
                                                                    'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(34, 48):
        if dataframe[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'post-seed', x, 0]], columns=['arm', 'type',
                                                                         'from_start',
                                                                         'pos'])
            dataframe = pd.concat([dataframe, new_row])
    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] > 47), 'from_start'] = 47

    for x in range(48, 53):
        if dataframe[(dataframe['type'] == 'loop') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['loop', 'loop', x, 0]], columns=['arm', 'type',
                                                                      'from_start',
                                                                      'pos'])
            dataframe = pd.concat([dataframe, new_row])

    dataframe = dataframe.groupby(['arm', 'type',
                                   'from_start'],
                                  as_index=False)[['pos']].sum()[['arm', 'type',
                                                                  'from_start',
                                                                  'pos']]
    return dataframe


def prepare_data_3p(df_temp):
    dataframe = df_temp[(df_temp['arm'] == '3p')].groupby(['arm', 'type', 'from_start'],
                                                          as_index=False)[['pos']].count()[['arm', 'type',
                                                                                            'from_start',
                                                                                            'pos']]

    try:

        dataframe_loop = df_temp[(df_temp['arm'] == 'loop') &
                                 (df_temp['from end'] > -6)].groupby(['arm', 'type', 'from end'],
                                                                     as_index=False)[['pos']].count()[['arm', 'type',
                                                                                                       'from end',
                                                                                                       'pos']]
        dataframe_loop['from_start'] = dataframe_loop['from end'].apply(
            lambda start: start + 6
        )
        dataframe_loop.drop('from end', inplace=True, axis=1)
    except KeyError:

        dataframe_loop = df_temp[(df_temp['arm'] == 'loop') &
                                 (df_temp['from_end'] > -6)].groupby(['arm', 'type', 'from_end'],
                                                                     as_index=False)[['pos']].count()[['arm', 'type',
                                                                                                       'from_end',
                                                                                                       'pos']]
        dataframe_loop['from_start'] = dataframe_loop['from_end'].apply(
            lambda start: start + 6
        )
        dataframe_loop.drop('from_end', inplace=True, axis=1)

    dataframe = pd.concat([dataframe, dataframe_loop], sort=False)

    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '3p'), 'from_start'] + 13

    dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '3p'), 'from_start'] + 6

    dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '3p'), 'from_start'] + 5

    dataframe.loc[(dataframe['type'] == 'flanking-3') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'flanking-3') & (dataframe['arm'] == '3p'), 'from_start'] + 27

    for x in range(28, 53):
        if dataframe[(dataframe['type'] == 'flanking-3') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'flanking-3', x, 0]], columns=['arm', 'type',
                                                                          'from_start',
                                                                          'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(14, 28):
        if dataframe[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'post-seed', x, 0]], columns=['arm', 'type',
                                                                         'from_start',
                                                                         'pos'])
            dataframe = pd.concat([dataframe, new_row])

    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] > 27), 'from_start'] = 27

    for x in range(7, 14):
        if dataframe[(dataframe['type'] == 'seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'seed', x, 0]], columns=['arm', 'type',
                                                                    'from_start',
                                                                    'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(6, 7):
        if dataframe[(dataframe['type'] == 'pre-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'pre-seed', x, 0]], columns=['arm', 'type',
                                                                        'from_start',
                                                                        'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(1, 6):
        if dataframe[(dataframe['type'] == 'loop') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['loop', 'loop', x, 0]], columns=['arm', 'type',
                                                                      'from_start',
                                                                      'pos'])
            dataframe = pd.concat([dataframe, new_row])

    dataframe['from_start'] = dataframe['from_start'].apply(lambda start: start * -1)

    dataframe = dataframe.groupby(['arm', 'type',
                                   'from_start'],
                                  as_index=False)[['pos']].sum()[['arm', 'type',
                                                                  'from_start',
                                                                  'pos']]

    return dataframe


@click.command()
@click.option('--output_folder', '-o')
def main(output_folder='/Users/martynaurbanek/Documents/private/repos/files/output/'):
    if output_folder is None:
        output_folder = '/Users/martynaurbanek/Documents/private/repos/files/output/'

    files = [[(x[0] + '/' + y, x[0].split('/')[-1].replace('DATA_', '')) for y in x[2]] for x in os.walk(output_folder)]
    flat_files = [file for sublist in files for file in sublist]
    csv_files = [file for file in flat_files if file[0].endswith('csv')]

    folders = [
            'ACC_22',
            'BLCA_22',
            'BRCA_22',
            'CESC_22',
            'CHOL_22',
            'COAD_22',
            'COAD2_22',
            'DLBC_22',
            'ESCA_22',
            'GBM_22',
            'KICH_22',
            'KIRC_22',
            'KIRP_22',
            'LAML_22',
            'LGG_22',
            'LIHC_22',
            'LUAD_22',
            'LUSC_new',
            'MESO_22',
            'OV_22',
            'PAAD_22',
            'PCPG1_22',
            'PCPG2_22',
            'PRAD_22',
            'READ1_222',
            'READ2_22',
            'SARC1_22',
            'SARC2_22',
            'SARC3_22',
            'SKCM_22',
            'STAD_22',
            'TCGT_22',
            'THCA_22',
            'THYM1_22',
            'THYM2_22',
            'UCEC_22',
            'UCS_22',
            'UVM_22',
            'pancancer'
        ]

    for all_mutations_filtered_mut_type_gene, folder_name in [
        file for file in csv_files if file[0].endswith('all_mutations.csv')
        and file[1] in folders
    ]:
        print(folder_name)

        df_temp = pd.read_csv(all_mutations_filtered_mut_type_gene)
        df_temp = df_temp[df_temp['mutation_type'] == 'subst']

        mutations = df_temp.shape[0]
        genes = df_temp['pre_name'].nunique()

        create_plot(df_temp, output_folder + 'plots/plot_miRNA_{}.pdf'.format(folder_name),
                    mutations, genes)

        # 5' dominant miRNAs

        df_5_dominant = df_temp[df_temp['balance'] == '5p']

        mutations = df_5_dominant.shape[0]
        genes = df_5_dominant['pre_name'].nunique()

        create_plot(df_5_dominant, output_folder + 'plots/plot_5p_balance_miRNA_{}.pdf'.format(folder_name),
                    mutations, genes)

        # 3p dominant miRNAs

        df_3_dominant = df_temp[df_temp['balance'] == '3p']

        mutations = df_3_dominant.shape[0]
        genes = df_3_dominant['pre_name'].nunique()

        create_plot(df_3_dominant, output_folder + 'plots/plot_3p_balance_miRNA_{}.pdf'.format(folder_name),
                    mutations, genes)

        # balanced miRNAs

        df_no_dominant = df_temp[df_temp['balance'] == 'both']

        mutations = df_no_dominant.shape[0]
        genes = df_no_dominant['pre_name'].nunique()

        create_plot(df_no_dominant, output_folder + 'plots/plot_balanced_3p_miRNA_{}.pdf'.format(folder_name),
                    mutations, genes)


if __name__ == "__main__":
    main()
