# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:23:10 2019

@author: Suha Naser-Khdour
"""
import os
import csv
from ete3 import Tree
import pandas as pd
from random import randint, sample, choice
from tqdm import tqdm


def root_tree(tree_file):
#root the tree in a midpoint node
    t = Tree(tree_file)
    outgroup = t.get_midpoint_outgroup() #find the midpoint node
    t.set_outgroup(outgroup) #set the midpoint node as an outgroup
    t.write(format=1, outfile= tree_file) #write in Newick formate
    return

def Create_Hetero_Files(model_file, iqtree_file, tree_file, site_info_file, parameter_file, trees_file, param_file_list, temp_file):
#root the tree in a midpoint node
    with open(tree_file, 'r') as input_file, open(temp_file, 'w') as output_file:
        for line in input_file:
            output_file.writelines(line)
    root_tree(temp_file) 
    
#extract the number of nodes in a tree
    with open(temp_file, 'r') as input_file:
        for line in input_file:
            nodes = line.count(':')

#extract model name and parameters from the best_model file
    for line in open(model_file):
        if 'F81' in line or 'JC' in line or 'JC69' in line:
            parameters = ['1', '1', '1', '1', '1', '1']
            if '+' in line:
                eind = line.index('+')
            elif ':' in line:
                eind = line.index(':')
            model = line[4:eind]
        elif '{' in line: 
            sind = line.index('{')
            eind = line.index('}')
            parameters = line[sind+1:eind].split(',')
            model = line[4:sind]
            if model == 'HKY' or model == 'HKY85' or model == 'K80' or model == 'K2P':
                parameters.insert(0, '1')
                parameters.insert(2, '1')
                parameters.insert(3, '1')
                parameters.insert(4, parameters[1])
                parameters.insert(5, '1')
            elif model == 'TN' or model == 'TN93' or model == 'TNe':
                parameters.insert(0, '1')
                parameters.insert(2, '1')
                parameters.insert(3, '1')
                parameters.insert(5, '1')
            elif model == 'K3P' or model == 'K81' or model == 'K81u':
                parameters.insert(0, '1')
                parameters.insert(3, parameters[2])
                parameters.insert(4, parameters[1])
                parameters.insert(5, '1')
            elif model == 'TVM' or model == 'TVMe':
                parameters.insert(4, parameters[1])
                parameters.insert(5, '1')
            elif model == 'TPM2' or model == 'TPM2u':
                parameters.insert(2, parameters[0])
                parameters.insert(3, '1')
                parameters.insert(4, parameters[1])
                parameters.insert(5, '1')
            elif model == 'TPM3' or model == 'TPM3u':
                parameters.insert(2, '1')
                parameters.insert(3, parameters[0])
                parameters.insert(4, parameters[1])
                parameters.insert(5, '1')
            elif model == 'TIM' or model == 'TIMe':
                parameters.insert(0, '1')
                parameters.insert(3, parameters[2])
                parameters.insert(5, '1')
            elif model == 'TIM2' or model == 'TIM2e':
                parameters.insert(2, parameters[0])
                parameters.insert(3, '1')
                parameters.insert(5, '1')
            elif model == 'TIM3' or model == 'TIM3e':
                parameters.insert(2, '1')
                parameters.insert(3, parameters[0])
                parameters.insert(5, '1')
            elif model == 'SYM' or model == 'GTR':
                parameters.insert(5, '1')
        if '+F{' in line:
            sind = line.index('+F{')
            line = line[sind+3:]
            eind = line.index('}')
            freq = line[:eind].split(',')
            is_FQ = False
        else:
            freq = ['0.25', '0.25', '0.25', '0.25']
            is_FQ = True
    parameters = parameters + freq
                
#extract proportions of number of taxa, variant sites and invariant sites from the iqtree file
    with open(iqtree_file, 'r') as input_file:
        reader = csv.reader(input_file, delimiter='\t')
        desired_rows = [row for row_number, row in enumerate(reader) if row_number == 34 or row_number == 37]
    df = pd.DataFrame(desired_rows)
    invariant = (round(int(df.iloc[1][5])/int(df.iloc[1][2]), 5))
    variant = (round(1 - invariant, 5))
    invariant = [str(invariant)] + freq
    variant = [str(variant)] + freq
    invariant = '\t'.join(invariant)
    variant = '\t'.join(variant)
    sind = str(df.iloc[0][0]).find(':')
    eind = str(df.iloc[0][0]).find('taxa')
    taxa = int(str(df.iloc[0][0])[sind+2:eind-1])

#write the site information file
    with open(site_info_file, 'w') as output_file:
        output_file.writelines('# Format: [name of site category]	[variant/invariant]	[proportion]	[freq(A)]	[freq(C)]	[freq(G)]	[freq(T)]\n')
        output_file.writelines('Constant_site	invariant\t'+invariant+'\n')
        output_file.writelines('Category_1	variant\t'+variant)    

#write the SRH parameter file    
    with open(parameter_file, 'w') as output_file:
        output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
        for i in range(0, taxa):
            output_file.writelines('T'+str(i)+'\t')
            output_file.writelines('%s\t' % p  for p in parameters)
            output_file.writelines('\n')
        for i in range(0, nodes - taxa):
            output_file.writelines('N'+str(i)+'\t')
            output_file.writelines('%s\t' % p  for p in parameters)
            output_file.writelines('\n')

#add names to internal nodes in the tree file
    temp_file2 = os.path.join(DirName,'temp2.txt') #create temporary file to store tree info
    with open(temp_file, 'r') as input_file, open(temp_file2, 'w') as output_file:
        i = 0
        for line in input_file:
            line = line.replace('tax', 'T')
            for word in line.split('):'):
                if i < (nodes - taxa):
                    output_file.writelines(word+')N'+str(i)+':')
                    i += 1
                else:
                    output_file.writelines(word)
                    break
    with open(temp_file2, 'r') as input_file, open(temp_file, 'w') as output_file:
        for line in input_file:
            output_file.writelines(line)
    os.remove(temp_file2) #delete the temporary file

#write the trees file   
    with open(temp_file, 'r') as input_file, open(trees_file, 'w') as output_file:
        i = 0
        output_file.writelines('# Format: [name of variant site category]	[newick tree format]\n')
        output_file.writelines('Category_1\t')                       
        for line in input_file:
            output_file.writelines(line)

#write the parameter list file   
    with open(param_file_list, 'w') as output_file:
        output_file.writelines('# Format: [name of variant site category]	[parameter file name]\n')
        output_file.writelines('Category_1\t')
        output_file.writelines(os.path.dirname(param_file_list)+'/parameter_1.txt')
    
    return nodes-taxa, parameters, is_FQ

def Stationary_Heterogeneous_matrix(temp_file, parameter_file, i, parameters):
#write a non-SRH parameter file
    parameters_2 = parameters[:]
    parameters_2[:-4] = sample(parameters_2[:-4], 6)
    while True:
        t = Tree(temp_file, format=1)
        taxa = len(t.get_leaves())
        r_node = t.search_nodes(name="N"+str(randint(0,i-1)))[0] #find a random node for splitting the tree
        detached_node = r_node.detach() #split the tree
        if taxa > 5:
            if len(detached_node.get_leaves()) > 2 and len(t.get_leaves()) > 2:
                break
        elif taxa <=5:
            if len(detached_node.get_leaves()) != 0 and len(t.get_leaves()) != 0:
                break
    param_1 = []
    for n in t.iter_descendants():
        param_1.append(n.name)
    param_2 = [r_node.name]
    for n in detached_node.iter_descendants():
        param_2.append(n.name)
    with open(parameter_file, 'w') as output_file:
        output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
        for i in param_1:
            output_file.writelines(i+'\t')
            output_file.writelines('%s\t' % p  for p in parameters)
            output_file.writelines('\n')
        for i in param_2:
            output_file.writelines(i+'\t')
            output_file.writelines('%s\t' % p  for p in parameters_2)
            output_file.writelines('\n')
    return

def nonStationary_Heterogeneous_matrix(temp_file, parameter_file, i, parameters):
#write a non-SRH parameter file
    parameters[6:] = sample(parameters[6:], 4)
    parameters_2 = parameters[:]
    parameters_2[:-4] = sample(parameters_2[:-4], 6)
    while True:
        t = Tree(temp_file, format=1)
        taxa = len(t.get_leaves())
        r_node = t.search_nodes(name="N"+str(randint(0,i-1)))[0] #find a random node for splitting the tree
        detached_node = r_node.detach() #split the tree
        if taxa > 5:
            if len(detached_node.get_leaves()) > 2 and len(t.get_leaves()) > 2:
                break
        elif taxa <=5:
            if len(detached_node.get_leaves()) != 0 and len(t.get_leaves()) != 0:
                break
    param_1 = []
    for n in t.iter_descendants():
        param_1.append(n.name)
    param_2 = [r_node.name]
    for n in detached_node.iter_descendants():
        param_2.append(n.name)
    with open(parameter_file, 'w') as output_file:
        output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
        for i in param_1:
            output_file.writelines(i+'\t')
            output_file.writelines('%s\t' % p  for p in parameters)
            output_file.writelines('\n')
        for i in param_2:
            output_file.writelines(i+'\t')
            output_file.writelines('%s\t' % p  for p in parameters_2)
            output_file.writelines('\n')
    return

def nonStationary_Homogeneous_matrix(temp_file, parameter_file, i, parameters):
#write a non-SRH parameter file
    parameters[6:] = sample(parameters[6:], 4)
    t = Tree(temp_file, format=1)
    param_1 = []
    for n in t.iter_descendants():
        param_1.append(n.name)
    with open(parameter_file, 'w') as output_file:
        output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
        for i in param_1:
            output_file.writelines(i+'\t')
            output_file.writelines('%s\t' % p  for p in parameters)
            output_file.writelines('\n')
    return
        

if __name__ == '__main__': 
    rootDir = '/data/pandit/dna' #the root directory where all alignemnts are saved
    fns = [Stationary_Heterogeneous_matrix, nonStationary_Heterogeneous_matrix, nonStationary_Homogeneous_matrix]
    for DirName, subdirList, fileList in tqdm(os.walk(rootDir)):
        if 'data.nex.iqtree' in fileList:
            iqtree_file =  os.path.join(DirName,'data.nex.iqtree')
            tree_file =  os.path.join(DirName,'data.nex.treefile')
            model_file =  os.path.join(DirName,'data.nex.best_model.nex')
            site_info_file = os.path.join(DirName,'site_info_file.txt')
            parameter_file = os.path.join(DirName,'parameter_1.txt')
            trees_file = os.path.join(DirName,'trees.txt')
            param_file_list = os.path.join(DirName,'param_file_list.txt')
            temp_file = os.path.join(DirName,'temp.txt') #create temporary file to store tree info
            i, parameters, is_FQ, l = Create_Hetero_Files(model_file, iqtree_file, tree_file, site_info_file, parameter_file, trees_file, param_file_list,temp_file)
            if is_FQ:
                Stationary_Heterogeneous_matrix(temp_file, parameter_file, i, parameters)
            else:
                choice(fns)(temp_file, parameter_file, i, parameters)
            os.remove(temp_file) #delete the temporary file
            bashCommand = " ".join(["/home/suha/Hetero/Hetero2", trees_file, site_info_file, param_file_list, "-l", l])
            os.system(bashCommand) #run Hetero2 
