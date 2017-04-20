import numpy as np
import itertools as ite
from scipy.stats import chi2
import scipy as sp
from Bio.Nexus import Nexus
from Bio import AlignIO
import pandas as pd
from pathlib import Path
import time #only using time for timing/troubleshooting

#import scipy.stats as sp

def matrix(s1,s2):
    m= np.zeros((4,4))
    elog = [] #error log
    for i in range(len(s1)):
        #empty cases
        if s1[i] == "?":
            m = m
        elif s2[i] == "?":
            m = m
        elif s1[i] == "-":
            m = m
        elif s2[i] == "-":
            m = m
        #cases of A->A
        elif s1[i]==s2[i]:
            if s1[i]=="A":
                m[0,0]=m[0,0]+1
            elif s1[i]=="C":
                m[1,1]=m[1,1]+1
            elif s1[i]=="G":
                m[2,2]=m[2,2]+1
            elif s1[i]=="T":
                m[3,3]=m[3,3]+1
        #other cases
        elif s1[i]=="A":
            if s2[i]=="C":
                m[1,0]=m[1,0]+1
            elif s2[i]=="G":
                m[2,0]=m[2,0]+1
            elif s2[i]=="T":
                m[3,0]=m[3,0]+1
        elif s1[i]=="C":
            if s2[i]=="A":
                m[0,1]=m[0,1]+1
            elif s2[i]=="G":
                m[2,1]=m[2,1]+1
            elif s2[i]=="T":
                m[3,1]=m[3,1]+1
        elif s1[i]=="G":
            if s2[i]=="A":
                m[0,2]=m[0,2]+1
            elif s2[i]=="C":
                m[1,2]=m[1,2]+1
            elif s2[i]=="T":
                m[3,2]=m[3,2]+1
        elif s1[i]=="T":
            if s2[i]=="C":
                m[1,3]=m[1,3]+1
            elif s2[i]=="G":
                m[2,3]=m[2,3]+1
            elif s2[i]=="A":
                m[0,3]=m[0,3]+1           
        else:
            elog.append([s1[i],s2[i],"error in matrix"])
    return m, elog

def MPTS(m):
    """ inputs
            m: a 4x4 matrix of proportions
        outputs
            p: is a p-value for the matched pairs test of symmetry
    """
    s = 0.0
    for (i,j) in ite.product(range(0,4),range(0,4)):
        if i<j:
            n = (m[i,j]-m[j,i])**2
            d = m[i,j]+m[j,i]
            if float(d) != 0.:
                s = s+(float(n)/float(d)) 
                p = 1 - chi2.cdf(s,6.0)
            else:
                p = 'NA'

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
    for (i,j) in ite.product(range(0,3),range(0,3)):
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
    
def Test_aln(aln,dset):
    """inputs 
            charset_aln = alignment array of sites
        output
            p = array containing pvalues
    
    """
    
    aln_array = np.array([list(rec) for rec in aln], np.character)
    dat.charsets.keys() #these are the names to the CHARSETS in the .nex file, which you can iterate over in a for loop
    i = 0
    p = np.array(['Dataset','Charset','Test','Sp1','Sp2','p-value'],dtype='U14')
    for n in dat.charsets.keys():
        for q in ite.combinations(list(range(len(aln))),2): #iterating over all taxa for sites
            m, elog = matrix(aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())
            p=np.vstack([p,[dset,n,'MPTS',aln[q[0]].name,aln[q[1]].name,MPTS(m)]])
            p=np.vstack([p,[dset,n,'MPTMS',aln[q[0]].name,aln[q[1]].name,MPTMS(m)]])
        i = i+1
    return p

if __name__ == '__main__': 
    aln_path = input('input nex file here:')#'/Users/user/Documents/! ! 2017 ANU/Semester 1/SCNC2103/data reader/alignment.nex'
    start_time = time.time()
    dset=Path(aln_path).parts[-2]
    dat = Nexus.Nexus()

    dat.read(aln_path)
    
    aln = AlignIO.read(open(aln_path), "nexus")
    
    p = Test_aln(aln,dset)
    df = pd.DataFrame(p)
    df.to_csv("data.csv")
    print('process complete with no errors in', (time.time() - start_time))