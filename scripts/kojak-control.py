import pandas as pd
import sys

if len(sys.argv) <= 2:
    print 'Please specify input filename. First intra and then inter.'
    exit()

intra = pd.read_csv(sys.argv[1], sep='\t', index_col=False)
columns = intra.columns.tolist()
columns[-1] = 'Protein#1'
columns.append('Protein#2')

intra = pd.read_csv(sys.argv[1], sep='\t', index_col=False, skiprows=[0], header=None)
intra.columns = columns

inter = pd.read_csv(sys.argv[2], sep='\t', index_col=False, skiprows=[0], header=None)
inter.columns = columns

data = pd.concat([intra, inter])
data = data.sort_values(by=['scannr', 'Score'], ascending=False)
data = data.drop_duplicates(subset='scannr', keep='first')

intra = data[data['Protein#1'] == data['Protein#2']]
intra = intra.sort_values(by='Score', ascending=False)
intra = intra.reset_index()

inter = data[data['Protein#1'] != data['Protein#2']]
inter = inter.sort_values(by='Score', ascending=False)
inter = inter.reset_index()

#################
# Intra control
#################

intra['Decoy'] = intra['Protein#1'].str.startswith('DECOY_') & intra['Protein#2'].str.startswith('DECOY_')
intra['Decoy'] = intra['Decoy'].astype(int)
intra['CumDecoy'] = intra['Decoy'].cumsum()  # will count at the current position
intra['Total'] = pd.Series(range(1, len(intra)+1))

# control FDR using concatenate target-decoy approach
intra['FDR'] = intra['CumDecoy'] * 2 / intra['Total']
intra['FDR'] = intra['FDR'].clip(None, 1)

q_values = []
for i in range(0, len(intra)):
    q_values.append(intra[i:len(intra)]['FDR'].min())
intra['qValue'] = pd.Series(q_values)

intra = intra.drop(['index', 'Decoy', 'CumDecoy', 'Total', 'FDR'], axis=1)
filtered = intra[(~intra['Protein#1'].str.startswith('DECOY_')) & (~intra['Protein#2'].str.startswith('DECOY_'))]
filtered.to_csv(sys.argv[1] + '.controlled.csv', index=False)

#################
# Inter control
#################

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
    if inter['CumU'][i] < inter['CumF'][i]:
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
filtered = inter[(~inter['Protein#1'].str.startswith('DECOY_')) & (~inter['Protein#2'].str.startswith('DECOY_'))]
filtered.to_csv(sys.argv[2] + '.controlled.csv', index=False)
