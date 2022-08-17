# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:52:12 2019

@author: Suha Naser-Khdour
"""
import os
import sys
import numpy as np
from tqdm import tqdm
from Bio.Nexus import Nexus
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def partitionAln():
    if len(sys.argv) == 1:
        raise SystemExit('Error: Please provide the alignment path')
    if len(sys.argv) == 2:
        if sys.argv[1] == '-h':
            print('\nSyntax:\n')
            print('partitionAln.py <ALN_PATH>\n')
            print('<ALN_PATH> \t The path to the alignment.\n.')
            print('Notes:\n\t *The alignment should be in nexus format.\n\t *The alignment should contain a charset block with all the partitions.')
            return
        else:
            if os.path.isfile(sys.argv[1]):
                part(sys.argv[1])
            else:
                raise SystemExit('Error: No such file {}'.format(sys.argv[1]))
    elif len(sys.argv) > 2:
        raise SystemExit('Error: too many arguments')


def part(aln_path):
    '''
    Parameters
    ----------
    aln_path: path
         path to alignment file with Charset block in nexus format.
    Returns
    ----------
        all partition alignments in the same directory as the original alignment.
    '''
    aln = AlignIO.read(open(aln_path), "nexus")
    dat = Nexus.Nexus()
    dat.read(aln_path)
    aln_matrix = np.array([list(rec) for rec in aln])
    taxa_list = dat.get_original_taxon_order()
    for n in tqdm(dat.charsets.keys()):
        new_matrix = aln_matrix[:,dat.charsets[n]]
        new_aln = []
        for i in range(0, len(taxa_list)):
            new_aln.append(SeqRecord(Seq("".join(new_matrix[i,:])), id=taxa_list[i], annotations={"molecule_type": "DNA"}))
        AlignIO.write(MultipleSeqAlignment( new_aln[i] for i in range(0, len(taxa_list))), os.path.join(os.path.dirname(aln_path),n+'.nex'), "nexus")
    return

if __name__ == '__main__': 
    partitionAln()