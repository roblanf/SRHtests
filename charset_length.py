# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:23:10 2017

@author: Suha nasser
"""
import pandas as pd
import glob
import os
from Bio.Nexus import Nexus

table = pd.DataFrame(columns=['dataset', 'charset', 'length'])

'''
Calculate the length of the character set from the start position to the last position
'''
for f in glob.glob('/data/srh/SRHtests/datasets/*/alignment.nex'):

    dat = Nexus.Nexus()
    dat.read(f)

    for charset_name in dat.charsets.keys():                                 
        length = len(dat.charsets[charset_name])
        table = table.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(f)))), 'charset': charset_name, 'length': length}, ignore_index=True) 

table.to_csv('/data/srh/processed_data/SRH_tables/charset_length.csv')