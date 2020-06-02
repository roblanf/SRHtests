# -*- coding: utf-8 -*-
"""
Authors Wenqi Zhang, Suha Naser-Khdour, Rob Lanfear 
"""
import numpy as np
import itertools as ite
from scipy.stats import chi2
from Bio.Nexus import Nexus
from Bio import AlignIO
import pandas as pd
import math
import os
import shutil
from tqdm import tqdm

def nCr(n,r):
    '''
    The factorial function
    '''
    f = math.factorial
    return f(n) // f(r) // f(n-r)

def simMtx(a, x, y):
    '''
    divergence matrix
    '''
    a = np.array(list(a))
    x = np.array(list(x))
    y = np.array(list(y))
    ax = (x[:, None] == a[None, :]).astype(int)
    ay = (y[:, None] == a[None, :]).astype(int)
    return np.dot(ay.T, ax)

def MPTS(m):
    '''
    MaxSym test
    '''
    denominator = m+m.T
    off_diag_indices=np.triu_indices(len(denominator),1)
    numerator = np.power(m-m.T,2)
    numerator = np.squeeze(np.asarray(numerator[off_diag_indices]))
    denominator = np.squeeze(np.asarray(denominator[off_diag_indices]))
    nonzeros = np.where(denominator != 0)[0]
    if denominator[nonzeros].size != 0:
        s = np.sum(np.divide(numerator[nonzeros], denominator[nonzeros]))
        return float(s)
    else:
        return np.nan

def MPTS_df(m):
    '''
    degrees of freedom for the MaxSym test
    '''
    denominator = m+m.T
    off_diag_indices=np.triu_indices(len(denominator),1)
    denominator = np.squeeze(np.asarray(denominator[off_diag_indices]))
    i = np.count_nonzero(denominator)
    return int(i)

def MPTMS(m):
    '''
    MaxSym_mar test
    '''
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
    if np.linalg.matrix_rank(V) != V.shape[0]:
        return np.nan
    else:
        Vi=np.linalg.inv(V)
        s = (ut.dot(Vi)).dot(u)[0][0]
        return float(s)

def MPTIS(MPTSs,MPTSDF, MPTMSs):
    '''
    MaxSym_int test
    '''
    if isinstance(MPTSs,float) and isinstance(MPTMSs,float)==True:
        if (MPTSDF > 3):
            s = MPTSs-MPTMSs
        else:
            return np.nan
    else:
        return np.nan
    return float(s)

def pval(s,v):
    '''
    p-value of the chi-square test
    '''
    if math.isnan(s)==False and v>0:
        p = 1.-float(chi2.cdf(s,v))
        return p
    else:
        return np.nan

def Test_aln(aln,dset,dat):
    '''
    the matrix of all the pair-wise comparisons
    '''
    aln_array = np.array([list(rec) for rec in aln], np.character)
    dat.charsets.keys()
    i = 0
    no = nCr(len(aln),2)*3*len([len(v) for v in dat.charsets.keys()])+1
    p=np.empty([no,7],dtype='U22')
    p[0] = np.array(['dataset','Charset','Test','pvalue','d','Sp1','Sp2'])
    for n in tqdm(dat.charsets.keys()):
        for q in ite.combinations(list(range(len(aln))),2): #iterating over all taxa for sites
            m = simMtx('ACGT',aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())
            d = (np.sum(m)-sum(m[np.diag_indices(4)]))/np.sum(m)
            i = i+1
            p[i]=np.array([datas,n,'MPTS',pval(MPTS(m),MPTS_df(m)),d,aln[q[0]].name,aln[q[1]].name])
            i = i+1
            p[i]=np.array([datas,n,'MPTMS',pval(MPTMS(m),3),d,aln[q[0]].name,aln[q[1]].name])
            i = i+1
            p[i]=np.array([datas,n,'MPTIS',pval(MPTIS(MPTS(m),MPTS_df(m),MPTMS(m)),MPTS_df(m)-3),d,aln[q[0]].name,aln[q[1]].name])
    return p

def table(p):
    '''
    summary table of the pair-wise comparisons
    '''
    df = pd.DataFrame(p[1:], columns=p[0])
    df['pvalue'] = df['pvalue'].astype(float)
    df['d'] = df['d'].astype(int)
    df.drop(['Sp1','Sp2'], axis=1, inplace=True)
    df.dropna(inplace=True)
    T = df.groupby(['dataset','Charset','Test']).apply(lambda x: x.loc[x.d.idxmax(),['d','pvalue']]).reset_index()
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
    
    df = T.groupby(['Test']).get_group('MPTS')
    df = df.loc[df['pvalue']>=0.05]
    for i in df.Charset.unique():
        with open(MPTS_good_file,'a') as good_MPTS:
            good_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTS_all_file,'a') as all_MPTS:
            all_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTS')
    df = df.loc[df['pvalue']<0.05]
    for i in df.Charset.unique():
        with open(MPTS_bad_file,'a') as bad_MPTS:
            bad_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTS_all_file,'a') as all_MPTS:
            all_MPTS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTIS')
    df = df.loc[df['pvalue']>=0.05]
    for i in df.Charset.unique():
        with open(MPTIS_good_file,'a') as good_MPTIS:
            good_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTIS_all_file,'a') as all_MPTIS:
            all_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTIS')
    df = df.loc[df['pvalue']<0.05]
    for i in df.Charset.unique():
        with open(MPTIS_bad_file,'a') as bad_MPTIS:
            bad_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTIS_all_file,'a') as all_MPTIS:
            all_MPTIS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTMS')
    df = df.loc[df['pvalue']>=0.05]
    for i in df.Charset.unique():
        with open(MPTMS_good_file,'a') as good_MPTMS:
            good_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTMS_all_file,'a') as all_MPTMS:
            all_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    df = T.groupby(['Test']).get_group('MPTMS')
    df = df.loc[df['pvalue']<0.05]
    for i in df.Charset.unique():
        with open(MPTMS_bad_file,'a') as bad_MPTMS:
            bad_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
        with open(MPTMS_all_file,'a') as all_MPTMS:
            all_MPTMS.writelines(line for line in open(aln_path) if 'CHARSET '+i in line)
    
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
    aln_rootDir = '/data/srh/rawdata/'
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
                    dat = Nexus.Nexus()
                    dat.read(aln_path) 
                    aln = AlignIO.read(open(aln_path), "nexus")
                    p = Test_aln(aln,datas,dat)
                    TempDir = os.path.join(SRH_rootDir,datas,'Data')
                    if not os.path.exists(TempDir):
                        os.makedirs(TempDir)
                        os.chdir(TempDir)
                    df =pd.DataFrame(p[1:], columns=p[0])
                    T=table(p)
                    T.to_csv(os.path.join(TempDir,'MaxDiv.csv'))
                    df.to_csv(os.path.join(TempDir,'data2.csv'))
                    merged = pd.read_csv(os.path.join(TempDir,'merged.csv'))
                    merged = pd.merge(merged, T, how='outer', on=['dataset','Charset','Test'])
                    merged.drop(['isbad', 'Unnamed: 0'], axis=1, inplace=True)
                    merged.to_csv(os.path.join(TempDir,'merged2.csv'))
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
