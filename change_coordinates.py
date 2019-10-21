import pandas as pd

df = pd.read_csv('../files/files_for_mirnaome/hg19_knowngenes_ccds_oncodrivefml_mutacje_plik_corr',
                 header=None, sep='\t')
df[0] = df[0].apply(lambda x: 'chr' + str(x))
df = df[[0,1,2,4]]
print(df.head())

df.to_csv('../files/files_for_mirnaome/coordinates_ccds.bed', index=False, header=False,
          sep='\t')
