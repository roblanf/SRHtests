# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:00:09 2017

@author: Suha nasser
"""

import pandas as pd
import glob
import os

def charsets_MPTS(tablebinom):
    '''
    Calculate the Percentage of Bad character sets and Not-Bad character sets for the MPTS test
    '''
    csvin = pd.read_csv(tablebinom)
    csvin = csvin[csvin.Test == 'MPTS']
    bad = csvin.loc[csvin['p_binomial']<0.05]
    not_bad = csvin.loc[csvin['p_binomial']>=0.05]
    test_bad = len(bad)/len(csvin)*100
    test_not_bad = len(not_bad)/len(csvin)*100
    return test_bad, test_not_bad

def charsets_MPTIS(tablebinom):
    '''
    Calculate the Percentage of Bad character sets and Not-Bad character sets for the MPTIS test
    '''
    csvin = pd.read_csv(tablebinom)
    csvin = csvin[csvin.Test == 'MPTIS']
    bad = csvin.loc[csvin['p_binomial']<0.05]
    not_bad = csvin.loc[csvin['p_binomial']>=0.05]
    test_bad = len(bad)/len(csvin)*100
    test_not_bad = len(not_bad)/len(csvin)*100
    return test_bad, test_not_bad

def charsets_MPTMS(tablebinom):
    '''
    Calculate the Percentage of Bad character sets and Not-Bad character sets for the MPTMS test
    '''
    csvin = pd.read_csv(tablebinom)
    csvin = csvin[csvin.Test == 'MPTMS']
    bad = csvin.loc[csvin['p_binomial']<0.05]
    not_bad = csvin.loc[csvin['p_binomial']>=0.05]
    test_bad = len(bad)/len(csvin)*100
    test_not_bad = len(not_bad)/len(csvin)*100
    return test_bad, test_not_bad


all_charsets = pd.DataFrame(columns= ['dataset', 'test', 'Bad_Charsets %', 'Not_bad_Charsets %'])

for f in glob.glob('/data/srh/processed_data/SRH_tables/*/Data/tablebinom.csv'):
    MPTS_bad, MPTS_not_bad = charsets_MPTS(f)
    MPTIS_bad, MPTIS_not_bad = charsets_MPTIS(f)
    MPTMS_bad, MPTMS_not_bad = charsets_MPTMS(f)
    all_charsets = all_charsets.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(f))), 'test': 'MPTS', 'Bad_Charsets %': MPTS_bad, 'Not_bad_Charsets %': MPTS_not_bad}, ignore_index=True)
    all_charsets = all_charsets.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(f))), 'test': 'MPTIS', 'Bad_Charsets %': MPTIS_bad, 'Not_bad_Charsets %': MPTIS_not_bad}, ignore_index=True)
    all_charsets = all_charsets.append({'dataset': os.path.basename(os.path.dirname(os.path.dirname(f))), 'test': 'MPTMS', 'Bad_Charsets %': MPTMS_bad, 'Not_bad_Charsets %': MPTMS_not_bad}, ignore_index=True)
all_charsets.to_csv('/data/srh/tables/charsets_percentage.csv')