# -*- coding: utf-8 -*-
"""
Authors Wenqi Zhang, Suha Naser-Khdour, Rob Lanfear 
"""
import numpy as np
import itertools as ite
from scipy.stats import chi2
from scipy.stats import binom_test
import scipy as sp
from Bio.Nexus import Nexus
from Bio import AlignIO
import pandas as pd
import math
from pathlib import Path
import os
import shutil
from tqdm import tqdm

def nCr(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)

def simMtx(a, x, y):
    a = np.array(list(a))
    x = np.array(list(x))
    y = np.array(list(y))
    ax = (x[:, None] == a[None, :]).astype(int)
    ay = (y[:, None] == a[None, :]).astype(int)
    return np.dot(ay.T, ax)

def MPTS(m):
    d=(m+m.T)
    off_diag_indices=np.triu_indices(len(d),1)
    if 0 in d[off_diag_indices]:
        return np.nan
    else:
        numerator=(m-m.T)**2
        denominator=m+m.T
        return np.sum(numerator[off_diag_indices]/denominator[off_diag_indices])

def MPTMS(m):
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
            V[i,j]=r[i]+c[i]-2*m[i][i] #d_{i*}+d{*i}+2d{ii}
        elif i!=j:
            V[i,j]=-(m[i,j]+m[j,i])
    if sp.linalg.det(V) == 0:
        s=np.nan
    else:
        Vi=np.linalg.inv(V)
        s = (ut.dot(Vi)).dot(u)[0][0]
    return s

def MPTIS(MPTSs,MPTMSs):
    if isinstance(MPTSs,float) and isinstance(MPTMSs,float)==True:
        s = MPTSs-MPTMSs
    else:
        s=np.nan
    return s

def pval(sval,v):
    '''
    Gets a test statistic and outputs a pvalue for a chi squarred test with degrees of freedom v
    '''
    if math.isnan(sval)==False:
        p = 1.-float(chi2.cdf(sval,v))
    else:
        p = np.nan
    return p

def Test_aln(aln,dset,dat):
    aln_array = np.array([list(rec) for rec in aln], np.character)
    dat.charsets.keys()
    p = pd.DataFrame(columns=['Dataset','Charset','Test','Sp1','Sp2','pvalue'])
    for n in tqdm(dat.charsets.keys()):
        for q in ite.combinations(list(range(len(aln))),2): #iterating over all taxa for sites
            m = simMtx('ACGT',aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())
            p = p.append({'Dataset':dset,'Charset':n, 'Test':'MPTS', 'Sp1':aln[q[0]].name, 'Sp2':aln[q[1]].name, 'pvalue':pval(MPTS(m),6)}, ignore_index=True)
            p = p.append({'Dataset':dset,'Charset':n, 'Test':'MPTMS', 'Sp1':aln[q[0]].name, 'Sp2':aln[q[1]].name, 'pvalue':pval(MPTMS(m),3)}, ignore_index=True)
            p = p.append({'Dataset':dset,'Charset':n, 'Test':'MPTIS', 'Sp1':aln[q[0]].name, 'Sp2':aln[q[1]].name, 'pvalue':pval(MPTIS(MPTS(m),MPTMS(m)),3)}, ignore_index=True)
    return p
    
def table(p):
    Tests={'MPTS','MPTIS','MPTMS'}
    T = pd.DataFrame(columns=['Charset','Test','p<0.05','p>=0.05','NA','p_binomial'])
    for n in dat.charsets.keys():
        dfx = p.groupby(['Charset']).get_group(n)
        for m in Tests:
            M = dfx.groupby(['Test']).get_group(m)
            T = T.append({'Charset':n, 'Test':m,  'p<0.05':M.loc[M['pvalue']<0.05].count()['pvalue'], 'p>=0.05':M.loc[M['pvalue']>=0.05].count()['pvalue'], 'NA':sum(pd.isnull(M['pvalue'])), 'p_binomial': binom_test(M.loc[M['pvalue']<0.05].count()['pvalue'], M.loc[M['pvalue']<0.05].count()['pvalue'] + M.loc[M['pvalue']>=0.05].count()['pvalue'], 0.05, 'greater')}, ignore_index = True)
    return T

def init_partition_files(partition_file):
    with open(partition_file,'w') as f:
        f.writelines("#nexus\n")
        f.writelines("begin sets;\n")
    return

def end_partition_files(partition_file):    
    with open(partition_file,'a') as f:
        f.writelines("end;")
    return 

def partition_files(T,aln_path):
    init_partition_files(MPTS_all_file)
    init_partition_files(MPTIS_all_file)
    init_partition_files(MPTMS_all_file)
    init_partition_files(MPTS_good_file)
    init_partition_files(MPTIS_good_file)
    init_partition_files(MPTMS_good_file)
    init_partition_files(MPTS_bad_file)
    init_partition_files(MPTIS_bad_file)
    init_partition_files(MPTMS_bad_file)

    mpts = T.groupby(['Test']).get_group('MPTS')
    good = mpts.loc[mpts['p_binomial'] >= 0.05]['Charset'].values
    bad = mpts.loc[mpts['p_binomial'] < 0.05]['Charset'].values
    if len(good) > 0:
        for char in good:
            with open(MPTS_good_file,'a') as good_MPTS:
                good_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
            with open(MPTS_all_file,'a') as all_MPTS:
                all_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
    if len(bad) > 0:
        for char in bad:
            with open(MPTS_bad_file,'a') as bad_MPTS:
                bad_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
            with open(MPTS_all_file,'a') as all_MPTS:
                all_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
    
    mptis = T.groupby(['Test']).get_group('MPTIS')
    good = mptis.loc[mptis['p_binomial'] >= 0.05]['Charset'].values
    bad = mptis.loc[mptis['p_binomial'] < 0.05]['Charset'].values
    if len(good) > 0:
        for char in good:
            with open(MPTIS_good_file,'a') as good_MPTIS:
                good_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
            with open(MPTIS_all_file,'a') as all_MPTIS:
                all_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
    if len(bad) > 0:
        for char in bad:
            with open(MPTIS_bad_file,'a') as bad_MPTIS:
                bad_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
            with open(MPTIS_all_file,'a') as all_MPTIS:
                all_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
    
    mptms = T.groupby(['Test']).get_group('MPTMS')
    good = mptms.loc[mptms['p_binomial'] >= 0.05]['Charset'].values
    bad = mptms.loc[mptms['p_binomial'] < 0.05]['Charset'].values
    if len(good) > 0:
        for char in good:
            with open(MPTMS_good_file,'a') as good_MPTMS:
                good_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
            with open(MPTMS_all_file,'a') as all_MPTMS:
                all_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
    if len(bad) > 0:
        for char in bad:
            with open(MPTMS_bad_file,'a') as bad_MPTMS:
                bad_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
            with open(MPTMS_all_file,'a') as all_MPTMS:
                all_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+ char in line)
    
    end_partition_files(MPTS_all_file)
    end_partition_files(MPTIS_all_file)
    end_partition_files(MPTMS_all_file)
    end_partition_files(MPTS_good_file)
    end_partition_files(MPTIS_good_file)
    end_partition_files(MPTMS_good_file)
    end_partition_files(MPTS_bad_file)
    end_partition_files(MPTIS_bad_file)
    end_partition_files(MPTMS_bad_file)
    
    return
    
if __name__ == '__main__': 
    aln_rootDir = '/data/srh/SRHtests/datasets/'
    SRH_rootDir = '/data/srh/processed_data/SRH_tables/'
    IQtree_rootDir = '/data/srh/processed_data/IQtree/'
    for DirName, subdirList, fileList in os.walk(aln_rootDir):
        for fname in fileList:
            if(fname=="alignment.nex"):
                print(DirName)
                head_DirName, datas = os.path.split(DirName)
                aln_path = os.path.join(DirName,fname)

                # if this file exists, we already did this dataset successfully
                if not os.path.exists(os.path.join(IQtree_rootDir,datas,"MPTMS","Bad")):
                    dset=Path(aln_path).parts[-2]
                    dat = Nexus.Nexus()
                    dat.read(aln_path) 
                    aln = AlignIO.read(open(aln_path), "nexus")
                    p = Test_aln(aln,dset,dat)
                    TempDir = os.path.join(SRH_rootDir,datas,'Data')
                    if not os.path.exists(TempDir):
                        os.makedirs(TempDir)
                        os.chdir(TempDir)
                    T=table(p)
                    T.to_csv('tablebinom.csv')
                    p.to_csv('data.csv')
                    all_MPTS_path = os.path.join(IQtree_rootDir,datas,"MPTS","All")
                    if not os.path.exists(all_MPTS_path):
                        os.makedirs(all_MPTS_path)
                    good_MPTS_path = os.path.join(IQtree_rootDir,datas,"MPTS","Not_Bad")
                    if not os.path.exists(good_MPTS_path):
                        os.makedirs(good_MPTS_path)
                    bad_MPTS_path = os.path.join(IQtree_rootDir,datas,"MPTS","Bad")
                    if not os.path.exists(bad_MPTS_path):
                        os.makedirs(bad_MPTS_path)
                    all_MPTIS_path = os.path.join(IQtree_rootDir,datas,"MPTIS","All")
                    if not os.path.exists(all_MPTIS_path):
                        os.makedirs(all_MPTIS_path)
                    good_MPTIS_path = os.path.join(IQtree_rootDir,datas,"MPTIS","Not_Bad")
                    if not os.path.exists(good_MPTIS_path):
                        os.makedirs(good_MPTIS_path)
                    bad_MPTIS_path = os.path.join(IQtree_rootDir,datas,"MPTIS","Bad")
                    if not os.path.exists(bad_MPTIS_path):
                        os.makedirs(bad_MPTIS_path)
                    all_MPTMS_path = os.path.join(IQtree_rootDir,datas,'MPTMS','All')
                    if not os.path.exists(all_MPTMS_path):
                        os.makedirs(all_MPTMS_path)
                    good_MPTMS_path = os.path.join(IQtree_rootDir,datas,"MPTMS","Not_Bad")
                    if not os.path.exists(good_MPTMS_path):
                        os.makedirs(good_MPTMS_path)
                    bad_MPTMS_path = os.path.join(IQtree_rootDir,datas,"MPTMS","Bad")
                    if not os.path.exists(bad_MPTMS_path):
                        os.makedirs(bad_MPTMS_path)
                    MPTS_all_file = os.path.join(all_MPTS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(all_MPTS_path,'alignment.nex'))
                    MPTIS_all_file = os.path.join(all_MPTIS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(all_MPTIS_path,'alignment.nex'))
                    MPTMS_all_file = os.path.join(all_MPTMS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(all_MPTMS_path,'alignment.nex'))
                    MPTS_good_file = os.path.join(good_MPTS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(good_MPTS_path,'alignment.nex'))
                    MPTIS_good_file = os.path.join(good_MPTIS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(good_MPTIS_path,'alignment.nex'))
                    MPTMS_good_file = os.path.join(good_MPTMS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(good_MPTMS_path,'alignment.nex'))
                    MPTS_bad_file = os.path.join(bad_MPTS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(bad_MPTS_path,'alignment.nex'))
                    MPTIS_bad_file = os.path.join(bad_MPTIS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(bad_MPTIS_path,'alignment.nex'))
                    MPTMS_bad_file = os.path.join(bad_MPTMS_path,'partition.nex')
                    shutil.copy2(aln_path,os.path.join(bad_MPTMS_path,'alignment.nex'))
                    partition_files(T,aln_path)