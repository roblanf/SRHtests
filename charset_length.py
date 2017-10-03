# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:23:10 2017

@author: Suha nasser
"""
import pandas as pd
import glob
import os

table1 = pd.DataFrame(columns=['dataset', 'Charset', 'Start position', 'End position', 'Charset length'])
table2 = pd.DataFrame(columns=['dataset', 'Charset', 'Start position', 'End position', 'Charset length'])

'''
Calculate the length of a character set from the start position until the last position 
'''
for f in glob.glob('/data/srh/processed_data/IQtree/*/MPTIS/All/partition.nex'):
    with open (f, 'r') as nf:
        for line in nf:
            if "=" in line:
                sind = line.index('=')
                mind = line.index('-')
                if "\\" in line:
                    eind = line.index('\\')
                    start_site = line[sind+2:mind]
                    end_site = line[mind+1:eind]
                elif ";" in line:
                    eind = line.index(';')
                    start_site = line[sind+2:mind]
                    end_site = line[mind+1:eind]
                length = int(end_site) - int(start_site) + 1
                charset_name = line[9:sind-1]
                table1 = table1.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(f)))), 'Charset': charset_name, 'Start position': int(start_site), 'End position': int(end_site), 'Charset length': length}, ignore_index=True) 

'''
Calculate the length of a character set from the start position until the last position for the McCormack_2013 data set only
'''
f = '/data/srh/processed_data/IQtree/McCormack_2013/MPTIS/All/partition.nex'
with open (f, 'r') as nf:
    for line in nf:
        if "=" in line:
            sind = line.index('=')
            mind = line.index('-')
            if "\\" in line:
                eind = line.index('\\')
                start_site = line[sind+2:mind]
                end_site = line[mind+1:eind]
            elif ";" in line:
                eind = line.index(';')
                start_site = line[sind+2:mind]
                end_site = line[mind+1:eind]
            length = int(end_site) - int(start_site) + 1
            charset_name = line[12:sind-1]
            table2 = table2.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(f)))), 'Charset': charset_name, 'Start position': int(start_site), 'End position': int(end_site), 'Charset length': length}, ignore_index=True) 

table2 = pd.concat([table1, table2]).drop_duplicates(subset=['dataset','Start position', 'End position', 'Charset length'], keep='last').reset_index(drop=True)
table2.sort_values(by=['dataset', 'Start position'], ascending=[True, True], inplace=True)
table2 = table2.reset_index(drop=True)
table2.to_csv('/data/srh/processed_data/SRH_tables/charset_length.csv')