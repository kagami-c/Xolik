"""Split the original result file into intra/inter FDR control."""
import sys
import pandas as pd

if len(sys.argv) <= 1:
    print('Please specify input filename.')
    exit()

data = pd.read_csv(sys.argv[1])
data = data.drop('qValue', axis=1)

data['RawProtein#1'] = data['Protein#1'].apply(lambda x: x[6:] if x.startswith('DECOY_') else x)
data['RawProtein#2'] = data['Protein#2'].apply(lambda x: x[6:] if x.startswith('DECOY_') else x)

# intra = data[data['Protein#1'] == data['Protein#2']]
# inter = data[data['Protein#1'] != data['Protein#2']]

intra = data[data['RawProtein#1'] == data['RawProtein#2']]
inter = data[data['RawProtein#1'] != data['RawProtein#2']]
intra = intra.drop(['RawProtein#1', 'RawProtein#2'], axis=1)
inter = inter.drop(['RawProtein#1', 'RawProtein#2'], axis=1)

intra = intra.reset_index()
inter = inter.reset_index()


# re-control intra result
# intra['Decoy'] = intra['Protein#1'].str.startswith('DECOY_') & intra['Protein#2'].str.startswith('DECOY_')
# intra['Decoy'] = intra['Decoy'].astype(int)
# intra['CumDecoy'] = intra['Decoy'].cumsum()  # will count at the current position
# intra['Total'] = pd.Series(range(1, len(intra)+1))

# control FDR using concatenate target-decoy approach
# intra['FDR'] = intra['CumDecoy'] * 2 / intra['Total']
# intra['FDR'] = intra['FDR'].clip(None, 1)

# q_values = []
# for i in range(0, len(intra)):
#     q_values.append(intra[i:len(intra)]['FDR'].min())
# intra['qValue'] = pd.Series(q_values)

# intra = intra.drop(['index', 'Decoy', 'CumDecoy', 'Total', 'FDR'], axis=1)
# intra.to_csv(sys.argv[1] + '.intra.csv', index=False)
# filtered = intra[(~intra['Protein#1'].str.startswith('DECOY_')) & (~intra['Protein#2'].str.startswith('DECOY_'))]
# filtered.to_csv(sys.argv[1] + '.intra.filtered.csv', index=False)

intra['U'] = ((intra['Protein#1'].str.startswith('DECOY_') & ~intra['Protein#2'].str.startswith('DECOY_'))
              | (~intra['Protein#1'].str.startswith('DECOY_') & intra['Protein#2'].str.startswith('DECOY_')))
intra['F'] = intra['Protein#1'].str.startswith('DECOY_') & intra['Protein#2'].str.startswith('DECOY_')
intra['U'] = intra['U'].astype(int)
intra['F'] = intra['F'].astype(int)
intra['CumU'] = intra['U'].cumsum()
intra['CumF'] = intra['F'].cumsum()
intra['Total'] = pd.Series(range(1, len(intra)+1))

# calculate FDR
fdrs = []
for i in range(0, len(intra)):
    target = intra['Total'][i] - intra['CumU'][i] - intra['CumF'][i]
    target = float(target)
    if intra['CumU'][i] <= intra['CumF'][i]:
        fdrs.append(intra['CumF'][i] / target)
    else:
        fdrs.append((intra['CumU'][i] - intra['CumF'][i])/target)
intra['FDR'] = pd.Series(fdrs)
intra['FDR'] = intra['FDR'].clip(None, 1)

# calculate q value
q_values = []
for i in range(0, len(intra)):
    q_values.append(intra[i:len(intra)]['FDR'].min())
intra['qValue'] = pd.Series(q_values)

intra = intra.drop(['index', 'U', 'F', 'CumU', 'CumF', 'Total', 'FDR'], axis=1)
intra.to_csv(sys.argv[1] + '.intra.csv', index=False)
filtered = intra[(~intra['Protein#1'].str.startswith('DECOY_')) & (~intra['Protein#2'].str.startswith('DECOY_'))]
filtered.to_csv(sys.argv[1] + '.intra.filtered.csv', index=False)


# re-control inter results
inter['U'] = ((inter['Protein#1'].str.startswith('DECOY_') & ~inter['Protein#2'].str.startswith('DECOY_'))
              | (~inter['Protein#1'].str.startswith('DECOY_') & inter['Protein#2'].str.startswith('DECOY_')))
inter['F'] = inter['Protein#1'].str.startswith('DECOY_') & inter['Protein#2'].str.startswith('DECOY_')
inter['U'] = inter['U'].astype(int)
inter['F'] = inter['F'].astype(int)
inter['CumU'] = inter['U'].cumsum()
inter['CumF'] = inter['F'].cumsum()
inter['Total'] = pd.Series(range(1, len(inter)+1))

# calculate FDR
fdrs = []
for i in range(0, len(inter)):
    target = inter['Total'][i] - inter['CumU'][i] - inter['CumF'][i]
    target = float(target)
    if inter['CumU'][i] <= inter['CumF'][i]:
        fdrs.append(inter['CumF'][i] / target)
    else:
        fdrs.append((inter['CumU'][i] - inter['CumF'][i])/target)
inter['FDR'] = pd.Series(fdrs)
inter['FDR'] = inter['FDR'].clip(None, 1)

# calculate q value
q_values = []
for i in range(0, len(inter)):
    q_values.append(inter[i:len(inter)]['FDR'].min())
inter['qValue'] = pd.Series(q_values)

inter = inter.drop(['index', 'U', 'F', 'CumU', 'CumF', 'Total', 'FDR'], axis=1)
inter.to_csv(sys.argv[1] + '.inter.csv', index=False)
filtered = inter[(~inter['Protein#1'].str.startswith('DECOY_')) & (~inter['Protein#2'].str.startswith('DECOY_'))]
filtered.to_csv(sys.argv[1] + '.inter.filtered.csv', index=False)