# Load Xolik result from a list of result files

import os
import pandas as pd
from collections import defaultdict


# build dataframes in unified structure
def build_dataframe(data, engine_name, attribute, filename):
    df = pd.DataFrame(data[['ScanNum', 'Peptide#1', 'LinkSite#1', 'Peptide#2', 'LinkSite#2',
                            'Protein#1', 'Protein#2', 'Score']])
    df['Attribute'] = attribute
    df.insert(0, 'File', filename)
    df.insert(2, 'Engine', engine_name)
    return df


# Xolik interface
def load_xolik_result(path, prefix):
    """Load results that matches the prefix requirement one by one.
    """
    results = [name for name in os.listdir(path) if name.startswith(prefix) and name.endswith('.filtered.csv')]
    frames = []
    for filename in results:
        is_intra = None
        if filename.endswith('.csv.intra.filtered.csv'):  # intra
            is_intra = True
        elif filename.endswith('.csv.inter.filtered.csv'):
            is_intra = False
        else:
            raise Exception('Wrong filename')
        fullname = path + '/' + filename
        f = pd.read_csv(fullname)
        f = f[f['qValue'] <= 0.05]
        f['LinkSite#1'] = f['LinkSite#1'] + 1
        f['LinkSite#2'] = f['LinkSite#2'] + 1
        f['Peptide#1'] = f['Peptide#1'].apply(lambda x: x.split('.')[1])
        f['Peptide#2'] = f['Peptide#2'].apply(lambda x: x.split('.')[1])
        attrib = 'intra' if is_intra else 'inter'
        #name = filename[:-23]
        name = filename.split('.')[0]
        frames.append(build_dataframe(f, 'Xolik', attrib, name))
    return frames


# ECL interface
def load_ecl_result(path, prefix):
    results = [name for name in os.listdir(path) if name.startswith(prefix) and name.endswith('.csv')]
    frames = []
    for filename in results:
        is_intra = None
        if filename.endswith('.intra.csv'):  # intra
            is_intra = True
        elif filename.endswith('.inter.csv'):
            is_intra = False
        else:
            raise Exception('Wrong filename')
        fullname = path + '/' + filename
        f = pd.read_csv(fullname)
        f['qValue'] = f['q_value']
        f['ScanNum'] = f['scan_num']
        f = f[f['qValue'] <= 0.05]
        f['LinkSite#1'] = f['site_1']
        f['LinkSite#2'] = f['site_2']
        f['Peptide#1'] = f['peptide_1']
        f['Peptide#2'] = f['peptide_2']
        f['Protein#1'] = f['protein_1']
        f['Protein#2'] = f['protein_2']
        f['Score'] = f['score']
        attrib = 'intra' if is_intra else 'inter'
        name = filename.split('.')[0]
        frames.append(build_dataframe(f, 'ECL', attrib, name))
    return frames


# Kojak interface
def split_kojak_peptide_perc(row):
    record = row['peptide']
    middle = record.split('.')[1]
    words = middle.split('--')
    pep1 = words[0][:-1]
    pep2 = words[1][:-1]
    pep1_words = pep1.split('(')
    seq1 = pep1_words[0]
    site1 = int(pep1_words[1])
    pep2_words = pep2.split('(')
    seq2 = pep2_words[0]
    site2 = int(pep2_words[1])
    row['Peptide#1'] = seq1
    row['LinkSite#1'] = site1
    row['Peptide#2'] = seq2
    row['LinkSite#2'] = site2
    return row


def load_kojak_result(path, prefix):
    """.intra.filtered.tsv.dedup.csv"""
    results = [name for name in os.listdir(path) if name.startswith(prefix) and name.endswith('.dedup.csv')]
    frames = []
    for filename in results:
        is_intra = None
        if 'intra' in filename:
            is_intra = True
        elif 'inter' in filename:
            is_intra = False
        else:
            raise Exception("Wrong filename")
        fullname = path + '/' + filename
        f = pd.read_csv(fullname)
        f['qValue'] = f['q-value']
        f = f[f['qValue'] <= 0.05]
        f = f.apply(split_kojak_peptide_perc, axis=1)
        f['Protein#1'] = f['proteinId#1']
        f['Protein#2'] = f['proteinId#2']
        f['Score'] = f['score']
        attrib = 'intra' if is_intra else 'inter'
        name = filename[:-29]
        if len(f) == 0:
            continue
        frames.append(build_dataframe(f, 'Kojak', attrib, name))
    return frames


# pLink interface
def load_plink_result(filename):
    data = pd.read_csv(filename, sep='\t', skiprows=[0])
    if len(data) <= 0:
        return []
    first_half = data[data.index % 2 == 0]
    first_half.columns = ['#', 'Spectrum', 'PeptideNum', 'UniquePepNum', 'Samples', 'Score', 'Condition',
                          'S', 'SS', 'SSS']
    first_half = first_half.reset_index()
    first_half = first_half.drop(['index', 'S', 'SS', 'SSS'], axis=1)
    second_half = data[data.index % 2 == 1]
    second_half = second_half.reset_index()
    second_half = second_half.drop(['index', 'Score'], axis=1)
    data = pd.concat([first_half, second_half], axis=1)

    data['PeptideNum'] = data['PeptideNum'].astype(int)
    if len(data[data['PeptideNum'] > 1]) > 0:
        raise Exception('Existing multiple match for one spectrum')
    data = data.drop(['PeptideNum', 'UniquePepNum', 'Samples', 'Condition', '#', '*', '#,Spec', 'SampleID'], axis=1)

    def extract_scan_num(spec_name):
        words = spec_name.split('.')
        return int(words[1])

    data['ScanNum'] = data['Spectrum'].apply(extract_scan_num)

    def extract_filename(spec_name):
        words = spec_name.split('.')
        return words[0]

    data['File'] = data['Spectrum'].apply(extract_filename)

    def split_link(pname):
        pname = pname.split(':')[0]
        pname = pname.split('/')[0]
        words = pname.split('-')
        fwords = words[0].split('(')
        alpha = fwords[0]
        alpha_site = int(fwords[1][:-1])
        swords = words[1].split('(')
        beta = swords[0]
        beta_site = int(swords[1][:-1])
        return alpha, alpha_site, beta, beta_site

    def get_alpha(pname):
        (alpha, alpha_site, beta, beta_site) = split_link(pname)
        return alpha

    def get_alpha_site(pname):
        (alpha, alpha_site, beta, beta_site) = split_link(pname)
        return alpha_site

    def get_beta(pname):
        (alpha, alpha_site, beta, beta_site) = split_link(pname)
        return beta

    def get_beta_site(pname):
        (alpha, alpha_site, beta, beta_site) = split_link(pname)
        return beta_site

    data['Peptide#1'] = data['Peptide'].apply(get_alpha)
    data['LinkSite#1'] = data['Peptide'].apply(get_alpha_site)
    data['Peptide#2'] = data['Peptide'].apply(get_beta)
    data['LinkSite#2'] = data['Peptide'].apply(get_beta_site)
    data['Protein#1'] = data['Proteins'].apply(get_alpha)
    data['Protein#2'] = data['Proteins'].apply(get_beta)

    data = data.drop(['Mod_Sites', 'Calc_M', 'Delta_M', 'ppm', 'Peptide', 'Proteins'], axis=1)

    def fill_attrib(row):
        if row['Protein#1'] == row['Protein#2']:
            row['Attribute'] = 'intra'
        else:
            row['Attribute'] = 'inter'
        return row

    data = data.apply(fill_attrib, axis=1)

    df = pd.DataFrame(data[['ScanNum', 'Peptide#1', 'LinkSite#1', 'Peptide#2', 'LinkSite#2',
                            'Protein#1', 'Protein#2', 'Score', 'Attribute', 'File']])
    df['Engine'] = 'pLink'

    return [df]


# ECL2 interface
def load_ecl2_result(path, prefix):
    results = [name for name in os.listdir(path) if name.startswith(prefix) and name.endswith('.csv')]
    frames = []
    for filename in results:
        is_intra = None
        if filename.endswith('.intra.target.csv'):  # intra
            is_intra = True
        elif filename.endswith('.inter.target.csv'):
            is_intra = False
        else:
            raise Exception('Wrong filename')
        fullname = path + '/' + filename
        f = pd.read_csv(fullname)
        f['qValue'] = f['q_value']
        f = f[f['qValue'] <= 0.05]
        f['ScanNum'] = f['scan_num']

        def fill_rows(row):
            word = row['protein'].split(';')
            row['Protein#1'] = word[0]
            row['Protein#2'] = word[1][1:]
            peptide = row['peptide'].split('-')
            row['Peptide#1'] = peptide[0][1:-1]
            row['LinkSite#1'] = int(peptide[1])
            row['Peptide#2'] = peptide[2][1:-1]
            row['LinkSite#2'] = int(peptide[3])
            return row

        f = f.apply(fill_rows, axis=1)
        f['Score'] = f['score']
        if len(f) == 0:
            continue
        attrib = 'intra' if is_intra else 'inter'
        name = filename.split('.')[0]
        frames.append(build_dataframe(f, 'ECL2', attrib, name))
    return frames


# pLink2 interface
def load_plink2_result(path, prefix):
    results = [name for name in os.listdir(path) if name.startswith(prefix) and name.endswith('.csv')]
    frames = []
    for filename in results:
        fullname = path + '/' + filename
        try:
            f = pd.read_csv(fullname, header=None, skiprows=[0])
        except pd.errors.EmptyDataError:
            continue

        f.columns = ["Order", "Title", "Deprecated", "Charge", "Precursor_Mass", "Peptide", "Peptide_Type",
                     "Linker", "Peptide_Mass", "Modifications", "Evalue", "SVM_Score", "Precursor_Mass_Error(Da)",
                     "Precursor_Mass_Error(ppm)", "Proteins", "Protein_Type", "FileID", "LabelID", "Alpha_Matched",
                     "Beta_Matched", "Alpha_Evalue", "Beta_Evalue"]

        def get_scan(s):
            return int(s.split('=')[-1][:-1])

        f['ScanNum'] = f['Deprecated'].apply(get_scan)
        f['Score'] = f['Evalue']

        def parse(row):
            row['Peptide#1'] = row['Peptide'].split('-')[0].split('(')[0]
            row['Peptide#2'] = row['Peptide'].split('-')[1].split('(')[0]
            row['LinkSite#1'] = int(row['Peptide'].split('-')[0].split('(')[1][:-1])
            row['LinkSite#2'] = int(row['Peptide'].split('-')[1].split('(')[1][:-1])
            row['Protein#1'] = row['Proteins'].split('-')[0].split('(')[0].strip()
            row['Protein#2'] = row['Proteins'].split('-')[1].split('(')[0].strip()
            return row

        f = f.apply(parse, axis=1)

        def extract_filename(spec_name):
            return spec_name.split('.')[0]

        f['File'] = f['Title'].apply(extract_filename)

        def fill_attrib(row):
            if row['Protein#1'] == row['Protein#2']:
                row['Attribute'] = 'intra'
            else:
                row['Attribute'] = 'inter'
            return row

        f = f.apply(fill_attrib, axis=1)
        f['Engine'] = 'pLink2'
        frames.append(pd.DataFrame(f[['File', 'ScanNum', 'Engine',
            'Peptide#1', 'LinkSite#1', 'Peptide#2', 'LinkSite#2',
            'Protein#1', 'Protein#2', 'Score', 'Attribute']]))
    return frames


# overall adjustment
def adjust_representation(row):
    seq1 = row['Peptide#1']
    seq2 = row['Peptide#2']
    site1 = row['LinkSite#1']
    site2 = row['LinkSite#2']
    pro1 = row['Protein#1']
    pro2 = row['Protein#2']
    seq1 = seq1.replace('I', 'L')
    seq2 = seq2.replace('I', 'L')
    if seq1 > seq2 or (seq1 == seq2 and site1 > site2):
        seq1, seq2 = seq2, seq1
        site1, site2 = site2, site1
        pro1, pro2 = pro2, pro1
    pro1 = '|'.join(pro1.split('|')[:2])
    pro2 = '|'.join(pro2.split('|')[:2])
    row['Peptide#1'] = seq1
    row['Peptide#2'] = seq2
    row['LinkSite#1'] = site1
    row['LinkSite#2'] = site2
    row['Protein#1'] = pro1
    row['Protein#2'] = pro2
    return row


# generate venn diagram entries
def venn_spectra_entries(data):
    histogram = defaultdict(int)
    for fileindex in data['File'].unique():
        file_subset = data[data['File'] == fileindex]
        for scan_num in file_subset['ScanNum'].unique():
            subset = file_subset[file_subset['ScanNum'] == scan_num]
            flags = sorted(list(subset['Engine'].unique()))
            key = '-'.join(flags)
            histogram[key] += 1
    return histogram


def venn_psm_entries(data):
    histogram = defaultdict(int)
    for fileindex in data['File'].unique():
        file_subset = data[data['File'] == fileindex]
        for scan_num in file_subset['ScanNum'].unique():
            scan_subset = file_subset[file_subset['ScanNum'] == scan_num]
            for p1 in scan_subset['Peptide#1'].unique():
                p1_subset = scan_subset[scan_subset['Peptide#1'] == p1]
                for p2 in p1_subset['Peptide#2'].unique():
                    p2_subset = p1_subset[p1_subset['Peptide#2'] == p2]
                    flags = sorted(list(p2_subset['Engine'].unique()))
                    key = '-'.join(flags)
                    histogram[key] += 1
    return histogram


def venn_xpsm_entries(data):
    histogram = defaultdict(int)
    for fileindex in data['File'].unique():
        file_subset = data[data['File'] == fileindex]
        for scan_num in file_subset['ScanNum'].unique():
            scan_subset = file_subset[file_subset['ScanNum'] == scan_num]
            for p1 in scan_subset['Peptide#1'].unique():
                p1_subset = scan_subset[scan_subset['Peptide#1'] == p1]
                for p2 in p1_subset['Peptide#2'].unique():
                    p2_subset = p1_subset[p1_subset['Peptide#2'] == p2]
                    for l1 in p2_subset['LinkSite#1'].unique():
                        l1_subset = p2_subset[p2_subset['LinkSite#1'] == l1]
                        for l2 in l1_subset['LinkSite#2'].unique():
                            l2_subset = l1_subset[l1_subset['LinkSite#2'] == l2]
                            flags = sorted(list(l2_subset['Engine'].unique()))
                            key = '-'.join(flags)
                            histogram[key] += 1
    return histogram


def venn_peptide_entries(data):
    histogram = defaultdict(int)
    for p1 in data['Peptide#1'].unique():
        p1_subset = data[data['Peptide#1'] == p1]
        for p2 in p1_subset['Peptide#2'].unique():
            p2_subset = p1_subset[p1_subset['Peptide#2'] == p2]
            flags = sorted(list(p2_subset['Engine'].unique()))
            key = '-'.join(flags)
            histogram[key] += 1
    return histogram


def venn_xpeptide_entries(data):
    histogram = defaultdict(int)
    for p1 in data['Peptide#1'].unique():
        p1_subset = data[data['Peptide#1'] == p1]
        for p2 in p1_subset['Peptide#2'].unique():
            p2_subset = p1_subset[p1_subset['Peptide#2'] == p2]
            for l1 in p2_subset['LinkSite#1'].unique():
                l1_subset = p2_subset[p2_subset['LinkSite#1'] == l1]
                for l2 in l1_subset['LinkSite#2'].unique():
                    l2_subset = l1_subset[l1_subset['LinkSite#2'] == l2]
                    flags = sorted(list(l2_subset['Engine'].unique()))
                    key = '-'.join(flags)
                    histogram[key] += 1
    return histogram

