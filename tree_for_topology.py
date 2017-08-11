# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:23:10 2017

@author: Suha nasser
"""
import os
import glob

def tree_for_topology(tree_file): #Combines three IQ trees into one tree file
    with open(tree_file,'w') as tf:
        if os.path.exists(All_tree_file):
            tf.writelines(line for line in open(All_tree_file))
        if os.path.exists(Bad_tree_file):
            tf.writelines(line for line in open(Bad_tree_file))
        if os.path.exists(Not_Bad_tree_file):
            tf.writelines(line for line in open(Not_Bad_tree_file))
    return

def is_non_zero_file(fpath): #Checks if a tree file is not empty
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0
        
for partition in glob.glob('/data/srh/processed_data/IQtree/*/*'):
    All_tree_file = os.path.join(partition, 'All/partition.nex.treefile')
    Bad_tree_file = os.path.join(partition, 'Bad/partition.nex.treefile')
    Not_Bad_tree_file = os.path.join(partition, 'Not_Bad/partition.nex.treefile')
    tree_file = os.path.join(partition, 'alignment.treefile')
    tree_for_topology(tree_file)
    if not is_non_zero_file(tree_file): os.remove(tree_file) #Deletes the tree file if empty