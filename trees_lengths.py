# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:41:23 2017

@author: Suha nasser
"""

import pandas as pd
import os

def trees_length(IQtree):
    '''
    Create a table of the trees lengths from IQtree output file
    '''
    with open(temp_file, 'w') as new_file:
        new_file.writelines(line[43:] for line in open(IQtree) if 'Total tree length' in line)
    df = pd.read_csv(temp_file, names=['tree length'], index_col=False)
    df['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(IQtree))))
    df['test'] = os.path.basename(os.path.dirname(os.path.dirname(IQtree)))
    df['partition'] = os.path.basename(os.path.dirname(IQtree))
    return df

tree_rootDir = '/data/srh/processed_data/IQtree'
all_data = pd.DataFrame()
for DirName, subdirList, fileList in os.walk(tree_rootDir):
    if 'partition.nex.iqtree' in fileList:
        original_file =  os.path.join(DirName,'partition.nex.iqtree')
        temp_file = os.path.join(DirName,'topology.txt')
        all_data = all_data.append(trees_length(original_file), ignore_index=True)
        os.remove(temp_file)
all_data.to_csv('/data/srh/processed_data/tables/trees_lengths.csv', index=False)