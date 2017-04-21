#!/usr/bin/env python3

import numpy as np
import itertools 
from scipy.stats import chi2
import scipy as sp
from Bio.Nexus import Nexus
from Bio import AlignIO
import pandas as pd
from pathlib import Path
import time #only using time for timing/troubleshooting

def seq_matrix(alphabet="ACGT", seq1, seq2):
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

def MPTMS(m):
    """ inputs
            m: a 4x4 matrix of proportions
        outputs
            p: is a p-value for the matched pairs test of marginal symmetry
    """
    r = np.zeros((3))
    r[0]=np.sum(m[0])
    r[1]=np.sum(m[1])
    r[2]=np.sum(m[2])
    c = [sum(row[i] for row in m) for i in range(len(m[0]))]
    d = [r[0]-c[0],r[1]-c[1],r[2]-c[2]]
    ut = np.array([[d[0],d[1],d[2]]])
    u = ut.transpose()
    V = np.zeros((3,3))
    for (i,j) in itertools.product(range(0,3),range(0,3)):
        if i==j:
            V[i,j]=r[i]+c[i]+2*m[i][i] #d_{i*}+d{*i}+2d{ii}
        elif i!=j:
            V[i,j]=-(m[i,j]+m[j,i])
    if sp.linalg.det(V) == 0:
        p='NA'
    else:
        Vi=np.linalg.inv(V)
        s = (ut.dot(Vi)).dot(u)[0][0]
        #print(s)
        p = 1 - chi2.cdf(s,3.0)
    return p

def analyse_alignment(aln_path):
    aln_name = Path(aln_path).parts[-2]

    # read in the nexus alignment
    dat = Nexus.Nexus()
    dat.read(aln_path)
    
    # and again in another way
    aln = AlignIO.read(open(aln_path), "nexus")

    # turn the charset into an array    
    aln_array = np.array([list(rec) for rec in aln], np.character)

    p = np.array(['Dataset','Charset','Test','Sp1','Sp2','p-value'],dtype='U14')

    # iterate over all the charsets in the alignment
    for n in dat.charsets.keys():

        # iterate over all pairs of taxa, and calculate test statistics for each
        for q in itertools.combinations(list(range(len(aln))),2):

            # get the two seqs as a dot-product matrix
            m = seq_matrix(aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),
                           aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())

            p=np.vstack([p,
                        [aln_name,n,
                        'MPTS',
                        charset[q[0]].name,
                        charset[q[1]].name,
                        MPTS(m)]])

            #p=np.vstack([p,[aln_name,n,'MPTMS',charset[q[0]].name,charset[q[1]].name,MPTMS(m)]])

if __name__ == "__main__":
    folder = sys.argv[1]
    output_path = sys.argv[2]

    analyse_alignments(folder, output_path)
