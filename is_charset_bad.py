# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 11:01:57 2017

@author: Suha nasser
"""
import pandas as pd
import glob
import os
import numpy as np            

def bad_data(tablebinom):
    '''
    Find if the character set is rejecting the null hypothesis of symmetry
    '''
    csvin = pd.read_csv(tablebinom)
    csvin['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(tablebinom)))
    csvin['isbad'] = np.where(csvin.p_binomial<0.05, 1,0) #If the p-binomial of a character set is less than 0.05, then the character set isbad and it rejects the null hypothesis of symmetry
    csvin.to_csv(tablebinom, index=False)
    return

all_data = pd.DataFrame()

for f in glob.glob('/data/srh/processed_data/SRH_tables/*/Data/tablebinom.csv'):
    bad_data(f)
    df1 = pd.read_csv(f)
    all_data = all_data.append(df1,ignore_index=True)
all_data.to_csv('/data/srh/processed_data/SRH_tables/isbad.csv', index=False)
