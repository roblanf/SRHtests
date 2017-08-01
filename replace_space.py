# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:13:31 2017

@author: Suha nasser
"""
import glob
import shutil
import os

def replace_space(damaged_file): 
#fix alignment that contains incorrect spaces
    temp_file = '/data/raw_data/temp.nex'
    with open(temp_file,'w') as tf:
        with open(damaged_file,'r') as df:
            for line in df:
                if 'CHARSET	' in line:
                    line = line.replace('CHARSET	','CHARSET ')
                tf.writelines(line)
    shutil.copy2(temp_file, damaged_file)
    os.remove(temp_file)
    return


for f in glob.glob('/data/raw_data/*/alignment.nex'):
    replace_space(f)