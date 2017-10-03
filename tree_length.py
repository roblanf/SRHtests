# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:41:23 2017

@author: Suha nasser
"""

import pandas as pd
import os

def trees_length(IQtree):
    '''
    Creats a CSV file with the tree length as has been calculated by IQtree
    '''
    with open(temp_file, 'w') as new_file:
        new_file.writelines(line[43:] for line in open(IQtree) if 'Total tree length' in line)
    df = pd.read_csv(temp_file, names=['tree length'], index_col=False)
    df['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(IQtree))))
    df['test'] = os.path.basename(os.path.dirname(os.path.dirname(IQtree)))
    df['partition'] = os.path.basename(os.path.dirname(IQtree))
    return df

tree_rootDir = '/data/srh/processed_data/IQtree'
dfa = pd.DataFrame()
for DirName, subdirList, fileList in os.walk(tree_rootDir):
    if 'partition.nex.iqtree' in fileList:
        original_file =  os.path.join(DirName,'partition.nex.iqtree')
        temp_file = os.path.join(DirName,'topology.txt')
        dfa = dfa.append(trees_length(original_file), ignore_index=True)
        os.remove(temp_file)
dfa.to_csv(os.path.join(tree_rootDir, "tree_length.csv"), index=False)