#!/usr/bin/env python3

import numpy as np
import itertools 
from scipy.stats import chi2
import scipy as sp
from Bio.Nexus import Nexus
from Bio import AlignIO
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import glob
import os
import sys


def seq_matrix(seq1, seq2, alphabet="ACGT"):
    '''Build the dot-product of two strings. Details here: http://stackoverflow.com/questions/43511674/calculating-a-similarity-difference-matrix-from-equal-length-strings-in-python/43512150#43512150
    '''
    alphabet = np.array(list(alphabet))
    seq1 = np.array(list(seq1))
    seq2 = np.array(list(seq2))
    aseq1 = (seq1[:, None] == alphabet[None, :]).astype(int)
    aseq2 = (seq2[:, None] == alphabet[None, :]).astype(int)
    return np.dot(aseq2.T, aseq1)

def MPTS(m):
    """ inputs
            m: an nxn matrix of dot product of two sequences
        outputs
            p: is a p-value for the matched pairs test of symmetry
    """
    d=(m+m.T)
    off_diag_indices=np.triu_indices(len(d),1)
    if 0 in d[off_diag_indices]:
        return 'NA'
    else:
        numerator=(m-m.T)**2
        denominator=m+m.T
        s = np.sum(numerator[off_diag_indices]/denominator[off_diag_indices])

    # degrees of freedom
    l = m.shape[0]
    df = (l * (l - 1) / 2)
    p = 1 - chi2.cdf(s,df)

    return p

def analyse_alignments(folder, output_path):
    # Analyse all alignments in a folder

    outf=open(output_path,'ab')

    p = np.array(['Dataset','Charset','Test','Sp1','Sp2','p-value'],dtype='U14')
    np.savetxt(outf, p.reshape(1, p.shape[0]), delimiter=',', fmt='%s')

    for aln_path in glob.iglob(os.path.join(folder, '**', '*.nex'), recursive=True):
        print(aln_path)
        analyse_alignment(aln_path, outf)

    outf.close()

def analyse_alignment(aln_path, outf):

    aln_name = Path(aln_path).parts[-2]

    # read in the nexus alignment
    dat = Nexus.Nexus()
    dat.read(aln_path)
    
    # and again in another way
    aln = AlignIO.read(open(aln_path), "nexus")

    # turn the charset into an array    
    aln_array = np.array([list(rec) for rec in aln], np.character)


    # iterate over all the charsets in the alignment
    for n in tqdm(dat.charsets.keys()):

        # iterate over all pairs of taxa, and calculate test statistics for each
        for q in tqdm(itertools.combinations(list(range(len(aln))),2)):

            # get the two seqs as a dot-product matrix
            m = seq_matrix(aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),
                           aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())

            p=np.array( [aln_name,n,
                        'MPTS',
                        aln[q[0]].name,
                        aln[q[1]].name,
                        MPTS(m)])
            
            np.savetxt(outf, p.reshape(1, p.shape[0]), delimiter=',', fmt='%s')

if __name__ == "__main__":
    folder = sys.argv[1]
    output_path = sys.argv[2]

    analyse_alignments(folder, output_path)
