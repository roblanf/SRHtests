# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 10:06:35 2017

@author: Suha nasser
"""

import pandas as pd

df0 = pd.read_csv('/data/srh/tables/summary_charsets.csv')
df1 = pd.read_csv('/data/srh/tables/topology_tests.csv')
df2 = pd.read_csv('/data/srh/tables/trees_lengths.csv')
df3 = pd.DataFrame()
df4 = pd.DataFrame()
df5 = pd.DataFrame()
df6 = pd.DataFrame()
df7 = pd.DataFrame()
df8 = pd.DataFrame()

df3 = pd.merge(df1, df2, how='left', on=['dataset', 'test', 'partition'])
del df3['pKH']
del df3['pWKH']
del df3['cELW']
del df3['pAU']
df4 = pd.merge(df3, df0, how='left', on=['dataset', 'test'])
df5 = df4.loc[df3['partition'] == 'All', ['dataset', 'partition', 'test', 'tree', 'tree length', 'logL', 'deltaL', 'bpRELL', 'pSH', 'pWSH', 'Bad Charset length', 'Not-bad Charset length']]
df5['No_sites'] = df5['Bad Charset length'] + df5['Not-bad Charset length']
df5['No_substitutions'] = df5['No_sites']*df5['tree length']
del df5['Bad Charset length']
del df5['Not-bad Charset length']
df6 = df4.loc[df3['partition'] == 'Bad', ['dataset', 'partition', 'test', 'tree', 'tree length', 'logL', 'deltaL', 'bpRELL', 'pSH', 'pWSH', 'Bad Charset length']]
df6['No_sites'] = df6['Bad Charset length']
df6['No_substitutions'] = df6['No_sites']*df6['tree length']
del df6['Bad Charset length']
df7 = df4.loc[df3['partition'] == 'Not_Bad', ['dataset', 'partition', 'test', 'tree', 'tree length', 'logL', 'deltaL', 'bpRELL', 'pSH', 'pWSH', 'Not-bad Charset length']]
df7['No_sites'] = df7['Not-bad Charset length']
df7['No_substitutions'] = df7['No_sites']*df7['tree length']
del df7['Not-bad Charset length']
df8 = df5.append([df6, df7], ignore_index=True)
df8.sort_values(by=['dataset', 'partition', 'tree'], inplace=True)

df8.to_csv('/data/srh/tables/summary_trees.csv', index=False)