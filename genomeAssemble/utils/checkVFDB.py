#!/usr/bin/env /home/p/anaconda/anaconda3_5.2.0/bin/python3.6
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 13:38:45 2018

@author: x
"""

import os
import subprocess
from Bio import SeqIO
import pandas as pd
import argparse

BLASTP = '/home/p/blast/bin/blastp'
VFDB_PR = '/mnt/d/linux/W/NCBI/VFDB/VFDB_blast201806/VFDB_Pr'
VFDB_INFO = '/mnt/d/linux/W/NCBI/VFDB/VFDB_info.txt'
VFDB_HEAD = '/mnt/d/linux/W/NCBI/VFDB/VFDB_blast201806/VFDB_Pr.head'
VFDB_VF = '/mnt/d/linux/W/NCBI/VFDB/VFDB_blast201806/VFDB_Pr.VF'
KO_info = '/mnt/d/linux/W/NCBI/KEGG/KO.txt'

def VFDB_analysis(file_pr, outfolder = '.', VFDB_score_min = 0.6):
    
    folder = outfolder
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    command = '{BLASTP} -query {file_pr} -db {VFDB_PR} -outfmt 6 -num_threads 24 -out result_VFDB.tab -max_target_seqs 100 -task blastp-fast '.format(BLASTP=BLASTP, outfolder=outfolder, file_pr=file_pr, VFDB_PR=VFDB_PR)
    subprocess.run(command, shell=True)
    dc_ref = SeqIO.to_dict(SeqIO.parse(open(VFDB_PR),'fasta'))
    dc_ref_len = {k:len(v) for k,v in dc_ref.items()}
    dc_test = SeqIO.to_dict(SeqIO.parse(open(file_pr),'fasta'))
    dc_test_len = {k:len(v) for k,v in dc_test.items()}
    df_align = pd.read_csv('result_VFDB.tab', sep='\t',header=None, names=["query", "subject","identity","matchLength","missMatch","gap","qstart","qend","sstart", "send", "evalue", "bitScore"])
    df_align_select = df_align[df_align.loc[:,'identity']>70].copy()
    df_align_select['query_len'] = df_align_select.loc[:,'query'].apply(lambda x:dc_test_len[x])
    df_align_select['subject_len'] = df_align_select.loc[:,'subject'].apply(lambda x:dc_ref_len[x])
    df_align_select['match_score'] = df_align_select.apply(lambda x:x['matchLength']**2/(x['query_len'] * x['subject_len']),axis=1)
    df_align_select['confidence_score'] = df_align_select.apply(lambda x:x['matchLength']/x['query_len'],axis=1)
    df_align_select2 = df_align_select[(df_align_select['match_score']>VFDB_score_min) & (df_align_select['confidence_score']>VFDB_score_min)]
    df_align_select3 = df_align_select2.drop_duplicates(subset='query')
    df_vfdb_head = pd.read_csv(VFDB_HEAD, sep='\t', index_col=0)
    df_vfdb_VF = pd.read_csv(VFDB_VF, sep='\t', index_col=0)
    df_vfdb_info = pd.read_csv(VFDB_INFO, sep='\t', index_col=0)
    df_result = df_align_select3.join(df_vfdb_head,on='subject',how='left')
    df_result = df_result.join(df_vfdb_VF,on='subject',how='left')
    df_result = df_result.join(df_vfdb_info, on='VF_id',how='left')
    df_result.to_excel('VFDB_annotation.xlsx',index=False)


parser = argparse.ArgumentParser(description='VFDB analysis of protein sequences')
parser.add_argument('file_pr',help = 'input protein file')
parser.add_argument('-o', '--outfolder', help='output folder for the result', default='.')
parser.add_argument('-m', '--VFDB_score_min', help='min alignment score to determine a VFDB match, default=0.6', default=0.6)
f = parser.parse_args()
VFDB_analysis(file_pr=f.file_pr, outfolder=f.outfolder, VFDB_score_min=f.VFDB_score_min)