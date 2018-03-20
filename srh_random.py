# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 11:16:17 2018

@author: Suha nasser
"""

import pandas as pd
import numpy as np
import math
import os
import shutil


def nCr(n,r):
    '''Factorial function by Mark Tolonen
    '''
    f = math.factorial
    return f(n) // f(r) // f(n-r)

def nCr_20(table):
    '''
    returns a list of the datasets that have at least 20 different combinations of charsets
    '''
    lis = {}
    df = pd.read_csv((table))
    df.rename(columns={'No of Bad Charset':'b', 'No of Not-bad Charset':'g'}, inplace=True)
    df = df.dropna(axis=0, how='any', subset =['b','g'])
    df0 = df.loc[df['test']=='MPTS', ['dataset','b','g']]
    df0['b'] = df0['b'].astype("int")
    df0['g'] = df0['g'].astype("int")
    a = df0.values
    for i in range(len(a)):
        if nCr(a[i][1]+a[i][2],a[i][1]) > 20:
            lis.update({a[i][0]:(a[i][1],a[i][2])})
    return (lis)

def random_part(path,k,n_k):
    '''
    Generates 20 random partition files and thier complement parition files for the MPTS test only
    '''
    for i in range(1,21):
        lines = open(path).readlines()
        lines = lines[2:-1]
        bad_random_path = os.path.join(random_MPTS_path,'Bad','Bad'+str(i))
        if not os.path.exists(bad_random_path):
            os.makedirs(bad_random_path)
        bad_random_file = os.path.join(bad_random_path, 'partition.nex')   
        with open(bad_random_file,'a') as b:
            b.writelines("#nexus\n")
            b.writelines("begin sets;\n")
            b.writelines(np.random.choice(lines,k, replace=False))
            b.writelines("end;")
        shutil.copy2(os.path.join(aln_path,datas,'MPTS','All','alignment.nex'),os.path.join(bad_random_path,'alignment.nex'))
        lines = open(bad_random_file).readlines()
        lines = lines[2:-1]
        good_random_path = os.path.join(random_MPTS_path,'Not_Bad','Not_Bad'+str(i))
        if not os.path.exists(good_random_path):
            os.makedirs(good_random_path)
        good_random_file = os.path.join(good_random_path, 'partition.nex')
        with open(good_random_file,'a') as g:
            g.writelines("#nexus\n")
            g.writelines("begin sets;\n")
            g.writelines(line for line in open(path) if line not in open(bad_random_file))
            g.writelines("end;")
        shutil.copy2(os.path.join(aln_path,datas,'MPTS','All','alignment.nex'),os.path.join(good_random_path,'alignment.nex'))
    return


aln_path = '/data/srh/processed_data/IQtree/'
datasets = nCr_20('/data/srh/tables/summary_charsets.csv')
for DirName, subdirList, fileList in os.walk(aln_path):
    head_DirName, datas = os.path.split(DirName)
    if datas in datasets.keys():
        random_MPTS_path = os.path.join(aln_path,datas,'MPTS','Random')
        if not os.path.exists(random_MPTS_path):
            os.makedirs(random_MPTS_path)
        k, n_k = datasets.get(datas)
        random_part(os.path.join(aln_path,datas,'MPTS','All','partition.nex'),k,n_k)