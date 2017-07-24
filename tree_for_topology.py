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
        
tree_rootDir = 'C:/Users/Suha nasser/Documents/Processed_data/IQtree'
for partition in glob.glob('C:/Users/Suha nasser/Documents/Processed_data/IQtree/*/*'):
    All_tree_file = os.path.join(partition, 'All/tree_alignment.treefile')
    Bad_tree_file = os.path.join(partition, 'Bad/tree_alignment.treefile')
    Not_Bad_tree_file = os.path.join(partition, 'Not_Bad/tree_alignment.treefile')
    tree_file = os.path.join(partition, 'tree.treefile')
    tree_for_topology(tree_file)
    if not is_non_zero_file(tree_file): os.remove(tree_file) #Deletes the tree file if empty