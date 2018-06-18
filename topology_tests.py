# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:23:10 2017

@author: Suha nasser
"""
import os
import re
import csv
import pandas as pd

def trees_for_topology(IQtree): 
    '''
    Creats csv file for the topology tests
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
    if os.stat(topology_file_temp).st_size > 0:
        with open(topology_file_temp, 'r') as inf, open(topology_file, 'w') as outf:
            writer = csv.writer(outf, dialect='excel')
            for line in inf:
                line = splitting_pattern.sub(',', line.strip())
                writer.writerow(line.split(','))
    else:
        return
    csvin = pd.read_csv(topology_file)
    csvin['partition'] = os.path.basename(os.path.dirname(IQtree))
    csvin['test'] = os.path.basename(os.path.dirname(os.path.dirname(IQtree)))
    csvin['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(IQtree))))
    csvin['Tree'].replace([1,2,3], ['All','pass','fail'], inplace=True) #Replaces the numbers with the relevant tree
    csvin.rename(columns={'Tree': 'tree'}, inplace=True)
    csvin['logL'] = -csvin['logL']
    return csvin

tree_rootDir = '/data/srh/processed_data/IQtree'
all_data = pd.DataFrame()
for DirName, subdirList, fileList in os.walk(tree_rootDir):
    if 'partition.nex.iqtree' in fileList:
        original_file =  os.path.join(DirName,'partition.nex.iqtree')
        topology_file_temp = os.path.join(DirName,'topology.txt')
        topology_file = os.path.join(DirName,'topology.csv')
        all_data = all_data.append(trees_for_topology(original_file), ignore_index=True)
        os.remove(topology_file_temp)
        if os.path.isfile(topology_file):
            if os.stat(topology_file).st_size == 0:
                os.remove(topology_file)
all_data.to_csv('/data/srh/tables/topology_tests.csv')