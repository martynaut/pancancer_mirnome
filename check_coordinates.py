import pandas as pd

df_1600 = pd.read_csv(
    '/Users/martynaurbanek/Documents/private/repos/files/file_for_mirnaome_test/1600_covered_mirnas_TCGA.txt',
    sep='\t')
df_1600 = df_1600[df_1600['HighCoverage'] == 1]
df_1600['name'] = df_1600['TargetID'].str.lower()

print(df_1600.head())

print(df_1600.shape)
coordinates = pd.read_csv(
'/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome_mirbase22.1/coordinates.bed',
    sep='\t', header=None
)
joined1 = df_1600.join(coordinates.set_index(3), on='name', how='outer', rsuffix='_mirbase')
print(joined1[joined1[1].isnull()])



coordinates2 = coordinates.join(df_1600.set_index('name'), on=3, how='left')
print(coordinates2[coordinates2['HighCoverage'] == 1].shape)
print(coordinates2.head())
print(coordinates2.shape)

df_old = pd.read_csv(
'/Users/martynaurbanek/Documents/private/repos/files/files_for_mirnaome/new_coordinates_all_02.bed',
    sep='\t', header=None
)

joined2 = coordinates.join(df_old.set_index(3), on=3, rsuffix='_old', how='left')
print(joined2[(joined2['1']==joined2['1_old']) & (joined2['2']==joined2['2_old'])].shape)
print(joined2[(joined2['1']!=joined2['1_old']) | (joined2['2']!=joined2['2_old'])].shape)

print(joined2[(joined2['1']!=joined2['1_old']) | (joined2['2']!=joined2['2_old'])].dropna().shape)
print(joined2[(joined2['1']!=joined2['1_old']) | (joined2['2']!=joined2['2_old'])].dropna())

