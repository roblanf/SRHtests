# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:38:28 2019

@author: Suha Naser-Khdour
"""

import re
import pandas as pd
import os
import csv
from ete3 import Tree
import sys

def parameters(IQtree):
    '''
    Creats a dataframe with all GTR parameters for each characterset
    This function works with iqtree files from IQ-TREE output
    the command line to estimate the GTR parameters from IQ-TREE is: iqtree -s alignmentFile -spp partitionFile -m GTR+I+G -n 0
    if the alignment file is not partitioned please use parameters2 function
    '''
    condition1 = False
    condition2 = False
    temp_file1 = os.path.join(os.path.dirname(IQtree), 'temp1.txt')
    temp_file2 = os.path.join(os.path.dirname(IQtree), 'temp2.txt')
    temp_df1 = os.path.join(os.path.dirname(IQtree), 'temp1.csv')
    temp_df2 = os.path.join(os.path.dirname(IQtree), 'temp2.csv')

    with open(temp_file1, 'w') as speed, open(temp_file2, 'w') as partitions:
        for line in open(IQtree, 'r'):
            if 'ID' and 'Model' and 'Speed' and 'Parameters' in line:
                condition1 = True
            if 'ID' and 'Name' and 'Type' and 'Seq' and 'Site' in line:
                condition2 = True

            if condition1 and not line.strip():
                condition1 = False
            if condition2 and not line.strip():
                condition2 = False

            if condition1:
                speed.writelines(line)
            if condition2:
                partitions.writelines(line)
            
    splitting_pattern = re.compile(r" {1,}")  #Splitting pattern if there is one or more spaces in a row
    if os.stat(temp_file1).st_size > 0:
        with open(temp_file1, 'r') as inf, open(temp_df1, 'w') as outf:
            writer = csv.writer(outf, dialect='excel')
            for line in inf:
                line = splitting_pattern.sub('\t', line.strip())
                writer.writerow(line.split('\t'))
    if os.stat(temp_file2).st_size > 0:
        with open(temp_file2, 'r') as inf, open(temp_df2, 'w') as outf:
            writer = csv.writer(outf, dialect='excel')
            for line in inf:
                line = splitting_pattern.sub('\t', line.strip())
                writer.writerow(line.split('\t'))
    else:
        return
    df = pd.read_csv(temp_df1)
    usecols  = [0,1,3,4]
    dftemp2 = pd.read_csv(temp_df2, usecols=usecols)
    df = pd.merge(df, dftemp2, how='outer', on=['ID'])
    df['dataset'] = os.path.basename(os.path.dirname(IQtree))
    df['P'] = df.Parameters.str.extract(r'\{(.*?)\}', expand=True)
    df['F'] = df.Parameters.str.extract(r'F\{(.*?)\}', expand=True)
    df['Invar'] = df.Parameters.str.extract(r'I\{(.*?)\}', expand=True)
    df['Gamma'] = df.Parameters.str.extract(r'G4\{(.*?)\}', expand=True)
    df[['F', 'C', 'G', 'T']] = df.F.str.split(',', expand=True)
    df[['P', 'A-G', 'A-T', 'C-G', 'C-T']] = df.P.str.split(',', expand=True)
    df['G-T'] = 1
    df.rename(columns={'Name': 'partition', 'F': 'A', 'P':'A-C', df.columns[5]:'sites', df.columns[6]:'taxa'}, inplace=True)
    df = df.drop([df.columns[0], df.columns[1], df.columns[3]], axis=1) 
    df = df[['dataset', 'partition', 'sites', 'taxa', 'A-C', 'A-G', 'A-T', 'C-G', 'C-T', 'G-T', 'A', 'C', 'G', 'T', 'Speed', 'Gamma', 'Invar']]

    os.remove(temp_file1)
    os.remove(temp_file2)
    os.remove(temp_df1)
    os.remove(temp_df2)
    return df


def parameters2(IQtree):
    '''
    Creats a dataframe with all GTR parameters for each characterset
    This function works with iqtree files from IQ-TREE output
    the command line to estimate the GTR parameters from IQ-TREE is: iqtree -s alignmentFile -m GTR+I+G -n 0
    if the alignment file is partitioned please use parameters function
    '''
    for line in open(IQtree, 'r'):
        if ('Input data:' in line) and ('sequences' in line) and ('sites' in line):
            seq = re.findall(r'Input data: ([\d.]+) sequences', line)
            site = re.findall(r'with ([\d.]+) nucleotide', line)
        if 'A-C:' in line:
            AC = re.findall(r'A-C: ([\d.]+)', line)
        if 'A-G:' in line:
            AG = re.findall(r'A-G: ([\d.]+)', line)
        if 'A-T:' in line:
               AT = re.findall(r'A-T: ([\d.]+)', line)
        if 'C-G:' in line:
            CG = re.findall(r'C-G: ([\d.]+)', line)
        if 'C-T:' in line:
            CT = re.findall(r'C-T: ([\d.]+)', line)
        if 'G-T:' in line:
            GT = re.findall(r'G-T: ([\d.]+)', line)
        if 'pi(A)' in line:
            A = re.findall(r'pi\(A\) = ([\d.]+)', line)
        if 'pi(C)' in line:
            C = re.findall(r'pi\(C\) = ([\d.]+)', line)
        if 'pi(G)' in line:
            G = re.findall(r'pi\(G\) = ([\d.]+)', line)
        if 'pi(T)' in line:
            T = re.findall(r'pi\(T\) = ([\d.]+)', line)
        if 'Proportion of invariable sites:' in line:
            invar = re.findall(r'sites: ([\d.]+)', line)
        if 'Gamma shape alpha:' in line:
            gamma = re.findall(r'alpha: *([\d.]+)', line)
            
    df = pd.DataFrame()
    df['sites'] = site
    df['taxa'] = seq
    df['invar'] = invar
    df['Gamma'] = gamma
    df['A'] = A
    df['C'] = C
    df['G'] = G
    df['T'] = T
    df['A-C'] = AC
    df['A-G'] = AG
    df['A-T'] = AT
    df['C-G'] = CG
    df['C-T'] = CT
    df['G-T'] = GT
    df['dataset'] = os.path.basename(os.path.dirname(IQtree))
    return df[['dataset', 'sites', 'taxa', 'A-C', 'A-G', 'A-T', 'C-G', 'C-T', 'G-T', 'A', 'C', 'G', 'T', 'Gamma', 'invar']]

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
    if (row.str.contains('16S').any() or row.str.contains('16s').any() or row.str.contains('28S').any() or row.str.contains('5p8S').any()
        or row.str.contains('18S').any() or row.str.contains('S7').any() or row.str.contains('lsu').any() or row.str.contains('LSU').any()):
        return 'rRNA'
    if (row.str.contains('Faircloth').any() or row.str.contains('McCormack').any() or row.str.contains('Moyle').any() 
        or row.str.contains('Prebus').any() or row.str.contains('Smith').any() or row.str.contains('Branstetter').any() 
        or row.str.contains('Crawford').any() or row.str.contains('Leache').any() or row.str.contains('uce').any()):
        return 'UCE'
    
    return 'other'

if __name__ == '__main__':
    print("In order to run this script all files must have the same name and extension and they should be saved in directories that have the datasets name. Please see an example below")
    diagram = Tree("((----->alignmentFileName.nex, ----->IQtreeFileName.iqtree)----->dataset1Dir, (----->alignmentFileName.nex, ----->IQtreeFileName.iqtree)----->dataset2Dir, (----->alignmentFileName.nex, ----->IQtreeFileName.iqtree)----->dataset3Dir)rootDir;", format=1)
    print(diagram.get_ascii(show_internal=True))
    proceed = input("do you want to proceed? Y/N\n")
    if proceed == 'Y':
        rootDir = '/data/Suha/GTR_parameters_dist/DNA/' #the rootDir name to the directories that contain the tree files
        IQtreeFileName = 'alignment.nex.iqtree' #the name of the iqtree file with .iqtree extension
        alignmentFileName = 'alignment.nex' #the name of the alignment file with extension
        parametersFile = 'GTRparam.csv' #the name of the GTR parameters output file with .csv extension
        
        df = pd.DataFrame()
    
        for DirName, subdirList, fileList in os.walk(rootDir):
            if IQtreeFileName in fileList:
                '''if you didn't allow different GTR models for each partition, please use parameters2 function instead of parametres function'''
                try:
                    iqtreeFile = os.path.join(DirName,IQtreeFileName)
                    df = df.append(parameters(iqtreeFile), ignore_index=True)
                except:
                    print('Could not add dataset: %s to the empirical distribution' % os.path.basename(DirName))

        charsetType = input("do you want to add a column for partition type: codon position, tRNA, rRNA, UCE or other? Y/N\n")
        if charsetType == 'Y':
                df['type'] = df.apply(nucleotide_type, axis=1)

        df.to_csv(os.path.join(rootDir, parametersFile))
    else:
        sys.exit()
