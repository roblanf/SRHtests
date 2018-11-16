# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 09:46:45 2018

@author: Suha Naser-Khdour
"""

from Bio import AlignIO
import os

def add_charset_block(input_file):
    '''
    input: alignment in nexus format
    output: alignment in nexus format with charset block
    '''
    aln = AlignIO.read(open(input_file), "nexus")
    l = len(aln[0])
    with open (input_file, 'a') as f:
        f.writelines("begin sets;\n")
        f.writelines("CHARSET    1st_pos = 1-"+str(l)+"\\3;\n")
        f.writelines("CHARSET    2nd_pos = 2-"+str(l)+"\\3;\n")
        f.writelines("CHARSET    3rd_pos = 3-"+str(l)+"\\3;\n")
        f.writelines("end;")

if __name__ == '__main__': 
    aln_rootDir = '/data/srh/SRHtests/datasets/' #Directory where all the alignment files are
    for DirName, subdirList, fileList in os.walk(aln_rootDir):
        for fname in fileList:
            if(fname=="alignment.nex"): #for convenience, all the alignments share the same name
                aln_path = os.path.join(DirName,fname)
                add_charset_block(aln_path)