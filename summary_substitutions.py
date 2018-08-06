# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 15:15:34 2018

@author: Suha Naser-Khdour
"""

import os
import re
import csv
import pandas as pd
from Bio.Nexus import Nexus
import glob

def genome(row):
    '''
    define charset's genome type: nuclear, mitochondrial, plastid or virus
    '''
    if row.str.contains('Broughton').any() and not row.str.contains('16s').any():
        return 'nuclear'
    if row.str.contains('Dornburg').any() and not row.str.contains('COI').any():
        return 'nuclear'
    if row.str.contains('Rightmyer').any() and row.str.contains('intron').any():
        return 'nuclear'
    if row.str.contains('Horn').any() and row.str.contains('EMB2765').any():
        return 'nuclear'
    if row.str.contains('Horn').any() and not row.str.contains('ITS').any() and not row.str.contains('nad1').any() and not row.str.contains('rps3').any() and not row.str.contains('EMB2765').any() and not row.str.contains('5p8S').any():
        return 'plastid'
    if row.str.contains('Sauquet').any() and not row.str.contains('its').any():
        return 'plastid'
    if row.str.contains('Oaks').any() and not row.str.contains('tRNA').any() and not row.str.contains('ND2').any() and not row.str.contains('ND3').any() and not row.str.contains('CYTB').any() and not row.str.contains('DLOOP').any():
        return 'nuclear'

    if row.str.contains('Cannon').any():
        return 'nuclear'
    if row.str.contains('Worobey').any():
        return 'virus'
    if row.str.contains('Wainwright').any():
        return 'nuclear'
    if row.str.contains('Fong').any():
        return 'nuclear'

    if row.str.contains('Faircloth').any():
        return 'nuclear'
    if row.str.contains('McCormack').any():
        return 'nuclear'
    if row.str.contains('Lartillot').any():
        return 'nuclear'
    if row.str.contains('Moyle').any():
        return 'nuclear'

    if row.str.contains('CO1').any() or row.str.contains('COI').any() or row.str.contains('co1').any():
        return 'mitochondrial'
    if row.str.contains('COII').any() or row.str.contains('CO2').any():
        return 'mitochondrial'
    if row.str.contains('CO3').any() or row.str.contains('co3').any():
        return 'mitochondrial'
    if row.str.contains('ND1').any() or row.str.contains('nd1').any():
        return 'mitochondrial'
    if row.str.contains('ND2').any() or row.str.contains('nd2').any():
        return 'mitochondrial'
    if row.str.contains('ND3').any() or row.str.contains('nd3').any():
        return 'mitochondrial'
    if row.str.contains('ND4').any() or row.str.contains('nd4').any():
        return 'mitochondrial'
    if row.str.contains('Cytb').any() or row.str.contains('CYTB').any() or row.str.contains('cytb').any():
        return 'mitochondrial'
    if row.str.contains('nad1').any():
        return 'mitochondrial'
    if row.str.contains('rps3').any():
        return 'mitochondrial'
    if row.str.contains('atp6').any():
        return 'mitochondrial'
    if row.str.contains('atp8').any():
        return 'mitochondrial'
    if row.str.contains('rps3').any():
        return 'mitochondrial'
    if row.str.contains('DLOOP').any():
        return 'mitochondrial'

    if row.str.contains('Phos').any():
        return 'nuclear'
    if row.str.contains('CAD').any() or row.str.contains('cad').any():
        return 'nuclear'
    if row.str.contains('ef1a').any() or row.str.contains('EF1a').any():
        return 'nuclear'
    if row.str.contains('Rag1').any() or row.str.contains('RAG1').any():
        return 'nuclear'
    if row.str.contains('Rag2').any():
        return 'nuclear'
    if row.str.contains('opsin').any():
        return 'nuclear'
    if row.str.contains('ACC').any():
        return 'nuclear'
    if row.str.contains('Cmos').any() or row.str.contains('cmos').any():
        return 'nuclear'
    if row.str.contains('Vim').any():
        return 'nuclear'
    if row.str.contains('PRLR').any():
        return 'nuclear'
    if row.str.contains('S7').any():
        return 'nuclear'
    if row.str.contains('H3').any():
        return 'nuclear'        
    if row.str.contains('ITS').any() or row.str.contains('its').any():
        return 'nuclear'        

    if row.str.contains('tRNA').any():
        return 'mitochondrial'
    if row.str.contains('16S').any() or row.str.contains('16s').any():
        return 'mitochondrial'
    if row.str.contains('18S').any():
        return 'nuclear'
    if row.str.contains('28S').any():
        return 'nuclear'
    if row.str.contains('5p8S').any():
        return 'nuclear'
    
    return 'other'

def nucleotide_type(row):
    '''
    define nucleotide's type: 1st codon position, 2nd codon position, 3rd codon position, tRNA, rRNA, UCE or other
    '''
    if row.str.contains('1st').any() or row.str.contains('pos1').any():
        return '1st_codon_position'
    if row.str.contains('2nd').any() or row.str.contains('2ns').any() or row.str.contains('pos2').any():
        return '2nd_codon_position'
    if row.str.contains('3rd').any() or row.str.contains('pos3').any():
        return '3rd_codon_position'
    if row.str.contains('tRNA').any():
        return 'tRNA'
    if row.str.contains('16S').any() or row.str.contains('18S').any() or row.str.contains('28S').any() or row.str.contains('5p8S').any():
        return 'rRNA'
    if row.str.contains('Faircloth').any() or row.str.contains('McCormack').any() or row.str.contains('Moyle').any():
        return 'UCE'
    
    return 'other'

def charsets_length(tree_rootDir):
    '''
    Create a dtaframe with the length of each characterset
    '''
    df = pd.DataFrame(columns=['dataset', 'charset', 'sites'])
    for f in glob.glob(os.path.join(tree_rootDir,'*','MPTS/All/alignment.nex')):
        dat = Nexus.Nexus()
        dat.read(f)
        for charset_name in dat.charsets.keys():
            sites = len(dat.charsets[charset_name])
            df = df.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(f)))), 'charset': charset_name, 'sites': sites}, ignore_index=True)
    
    return df

def speed(IQtree):
    '''
    Creats a dataframe with speeds for each characterset
    '''
    condition1 = False
    condition2 = False

    with open(temp_file1, 'w') as speed, open(temp_file2, 'w') as partitions:
        for line in open(IQtree, 'r'):
            if 'ID' and 'Model' and 'Speed' and 'Parameters' in line:
                condition1 = True
            if 'ID' and 'Name' and 'Type' and '#Seqs' and '#Sites' and '#Patterns' and '#Const_Sites' in line:
                condition2 = True

            if 'MAXIMUM LIKELIHOOD TREE' in line:
                condition1 = False
            if 'ModelFinder' in line:
                condition2 = False

            if condition1:
                line = line.replace(',', '') #Removes all ',' signs from the topology tests table
                speed.writelines(line)
            if condition2:
                partitions.writelines(line)

    splitting_pattern = re.compile(r" {1,}")  #Splitting pattern if there is one or more spaces in a row
    if os.stat(temp_file1).st_size > 0:
        with open(temp_file1, 'r') as inf, open(temp_df1, 'w') as outf:
            writer = csv.writer(outf, dialect='excel')
            for line in inf:
                line = splitting_pattern.sub(',', line.strip())
                writer.writerow(line.split(','))
    if os.stat(temp_file2).st_size > 0:
        with open(temp_file2, 'r') as inf, open(temp_df2, 'w') as outf:
            writer = csv.writer(outf, dialect='excel')
            for line in inf:
                line = splitting_pattern.sub(',', line.strip())
                writer.writerow(line.split(','))
    else:
        return

    df = pd.read_csv(temp_df1, usecols=['ID', 'Speed'])
    dftemp2 = pd.read_csv(temp_df2, usecols=['ID', 'Name'])
    df = pd.merge(df, dftemp2, how='left', on=['ID'])
    df['partition'] = os.path.basename(os.path.dirname(IQtree))
    df['partition'].replace(['Bad','Not_Bad'], ['fail','pass'], inplace=True) #Replaces the names of the subsets
    df['test'] = os.path.basename(os.path.dirname(os.path.dirname(IQtree)))
    df['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(IQtree))))

    return df

def taxa(IQtree):
    '''
    Create a dataframe of the trees lengths
    '''
    df = pd.DataFrame(columns=['dataset', 'test', 'partition', '#taxa'])
    for line in open(IQtree, 'r'):
        if 'taxa' in line:
            st = line.index(':') + 2
            end = line.index('taxa') - 1
            no_taxa = line[st:end]
    df = df = df.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(IQtree)))), 'test': os.path.basename(os.path.dirname(os.path.dirname(IQtree))), 'partition': os.path.basename(os.path.dirname(IQtree)), '#taxa':no_taxa}, ignore_index=True)
    df['partition'].replace(['Bad','Not_Bad'], ['fail','pass'], inplace=True) #Replaces the names of the subsets
    
    return df

def trees_length(IQtree):
    '''
    Create a dataframe of the trees lengths
    '''
    with open(temp_file1, 'w') as lengths:
        lengths.writelines(line[43:] for line in open(IQtree) if 'Total tree length' in line)
    df = pd.read_csv(temp_file1, names=['tree_length'], index_col=False)
    df['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(IQtree))))
    df['test'] = os.path.basename(os.path.dirname(os.path.dirname(IQtree)))
    df['partition'] = os.path.basename(os.path.dirname(IQtree))
    
    return df

if __name__ == '__main__': 
    tree_rootDir = '/data/srh/processed_data/IQtree'
    all_data = pd.DataFrame()
    tree_length = pd.DataFrame()
    No_taxa = pd.DataFrame()
    for DirName, subdirList, fileList in os.walk(tree_rootDir):
        if 'partition.nex.iqtree' in fileList:
            original_file =  os.path.join(DirName,'partition.nex.iqtree')
            temp_file1 = os.path.join(DirName,'speed.txt')
            temp_file2 = os.path.join(DirName,'partitions.txt')
            temp_df1 = os.path.join(DirName,'speed.csv')
            temp_df2 = os.path.join(DirName,'partitions.csv')
            all_data = all_data.append(speed(original_file), ignore_index=True)
            tree_length = tree_length.append(trees_length(original_file), ignore_index=True)
            No_taxa = No_taxa.append(taxa(original_file), ignore_index=True)
            os.remove(temp_file1)
            os.remove(temp_file2)
            os.remove(temp_df1)
            os.remove(temp_df2)
   
    tree_length['partition'].replace(['Bad','Not_Bad'], ['fail','pass'], inplace=True) #Replaces the names of the subsets
    all_data = pd.merge(all_data, tree_length, how='left', on=['dataset', 'test', 'partition'])
    all_data.rename(columns={'Name': 'charset'}, inplace=True)
    all_data = pd.merge(all_data, No_taxa, how='left', on=['dataset', 'test', 'partition'])
    charsets = charsets_length(tree_rootDir)
    all_data = pd.merge(all_data, charsets, how='left', on=['dataset', 'charset'])
    all_data['genome'] = all_data.apply(genome, axis=1)
    all_data['type'] = all_data.apply(nucleotide_type, axis=1)

    all_data.to_csv('/data/srh/tables/summary_substitutions.csv')