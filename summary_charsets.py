# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:12:55 2017

@author: Suha nasser
"""

import pandas as pd

df0 = pd.read_csv('/data/srh/tables/charsets_percentage.csv')
df1 = pd.read_csv('/data/srh/tables/is_charset_bad.csv')
df2 = pd.read_csv('/data/srh/tables/charset_length.csv')
df3 = pd.DataFrame()
df4 = pd.DataFrame()
df5 = pd.DataFrame()
df3 = pd.merge(df1, df2, how='left', left_on=['dataset', 'Charset'], right_on=['dataset', 'charset'])
df6 = df3.loc[df3['isbad'] == 0, ['dataset', 'Test', 'length']]
df6 = df6.groupby(['dataset', 'Test']).size().reset_index()
df7 = df3.loc[df3['isbad'] == 1, ['dataset', 'Test', 'length']]
df7 = df7.groupby(['dataset', 'Test']).size().reset_index()
df3 = df3.groupby(['dataset', 'Test', 'isbad']).sum().reset_index()
df4 = df3.loc[df3['isbad'] == 0, ['dataset', 'Test', 'length']]
df5 = df3.loc[df3['isbad'] == 1, ['dataset', 'Test', 'length']]
df4.rename(columns={'length': 'pass Charset length'}, inplace=True)
df5.rename(columns={'length': 'fail Charset length'}, inplace=True)
df4 = pd.merge(df4, df6, how = 'left', left_on=['dataset', 'Test'], right_on=['dataset', 'Test'])
df5 = pd.merge(df5, df7, how = 'left', left_on=['dataset', 'Test'], right_on=['dataset', 'Test'])
df0 = pd.merge(df0, df5, how = 'left', left_on=['dataset', 'test'], right_on=['dataset', 'Test'])
df0 = pd.merge(df0, df4, how = 'left', left_on=['dataset', 'test'], right_on=['dataset', 'Test'])
df0.rename(columns={'Test_x': 'Test_fail Charsets', 'Test_y': 'Test_pass Charsets', '0_y': 'No of pass Charset', '0_x': 'No of fail Charset'}, inplace=True)
df0['Avg(fail Charset length)'] = df0['fail Charset length']/df0['No of fail Charset']
df0['Avg(pass Charset length)'] = df0['pass Charset length']/df0['No of pass Charset']
df0['Delta length (fail-pass)'] = df0['Avg(fail Charset length)'] - df0['Avg(pass Charset length)']
df0.to_csv('/data/srh/tables/summary_charsets.csv', index=False)