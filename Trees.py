# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:23:10 2017

@author: Suha nasser
"""
import os
import re
import csv
import glob
import pandas as pd

def trees_for_topology(IQtree, tree_name, test_name, dataset_name):
    '''
    A function that creats a CSV table for the topology tests results
    '''
    condition = False
    with open(topology_file_temp, 'w') as new_file:
        for line in open(IQtree):
            if 'deltaL' in line: #Strats copying the IQtree file from the topology tests line
                condition = True
            if 'logL difference from the maximal logl in the set' in line: #Stops copying the IQtree file after the end of the topology tests lines
                condition = False
            if condition:
                line = line.replace('+', '') #Removes all '+' signs from the topology tests table
                line = line.replace('-', '') #Removes all '-' signs from the topology tests table
                new_file.writelines(line)
    splitting_pattern = re.compile(r" {1,}")  #Splitting pattern if there is one or more spaces in a row
    with open(topology_file_temp, 'r') as inf, open(topology_file, 'w') as outf:
        writer = csv.writer(outf, dialect='excel')
        for line in inf:
            line = splitting_pattern.sub(',', line.strip())
            writer.writerow(line.split(','))
    csvin = pd.read_csv(topology_file)
    csvin['Partition'] = tree_name
    csvin['Test'] = test_name
    csvin['Dataset'] = dataset_name
    csvin['Tree'].replace([1,2,3], ['All','Bad','Not_Bad'], inplace=True) #Replaces the numbers with the relevant tree
    csvin.to_csv(topology_file, index=False)
    return

tree_rootDir = '/data/srh/processed_data/IQtree/'
for DirName, subdirList, fileList in os.walk(tree_rootDir):
    if 'alignment.iqtree' in fileList:
        original_file =  os.path.join(DirName,'alignment.iqtree')
        topology_file_temp = os.path.join(DirName,'topology.txt')
        topology_file = os.path.join(DirName,'topology.csv')
        tree_name = os.path.basename(DirName)
        test_name = os.path.basename(os.path.dirname(DirName))
        dataset_name = os.path.basename(os.path.dirname(os.path.dirname(DirName)))
        trees_for_topology(original_file, tree_name, test_name, dataset_name)
        os.remove(topology_file_temp)
all_data = pd.DataFrame()
for f in glob.glob('/data/srh/processed_data/IQtree//*/*/*/topology.csv'): #Copies all the topology files into one file 
    df = pd.read_csv(f)
    all_data = all_data.append(df,ignore_index=True)
all_data.to_csv(os.path.join(tree_rootDir, "final_topology.csv"))