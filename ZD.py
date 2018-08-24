#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 17:22:58 2018

@author: x
"""

import glob
folder_abi = '/mnt/d/linux/W/moon/MoonZZD/abi/'
files = glob.glob(folder_abi + '**/*.ab1', recursive=True)
print(len([e for e in files if '.ab1' in e]))
print(len([e for e in files if '-27F.ab1' in e]))
print(len([e for e in files if '-1492R.ab1' in e]))


import pandas as pd
df = pd.read_csv('/mnt/d/linux/W/moon/MoonZZD/nameNotes.txt', sep='\t',dtype=str)

import os
folders = os.listdir(folder_abi)

folder_mnh = '/mnt/d/linux/W/moon/MoonZZD/abi_mnh/'
for _folder in folders:
    os.mkdir(os.path.join(folder_mnh,_folder))

df_withna = df[df.Nameid.isna() | df.folder.isna()]
df_withna.to_excel('/mnt/d/linux/W/moon/MoonZZD/nameNotes_withna.xlsx')
df_good = df.dropna()
df_good = df_good.copy()

df_good['file_ori_27F'] = df_good.apply(lambda x:os.path.join(folder_abi,x['folder'],x['Nameid']+'-27F.ab1'),axis=1)
df_good['file_ori_1492R'] = df_good.apply(lambda x:os.path.join(folder_abi,x['folder'],x['Nameid']+'-1492R.ab1'),axis=1)
df_good['file_des_27F'] = df_good.apply(lambda x:os.path.join(folder_mnh,x['folder'],x['Moonid']+'-27F.ab1'),axis=1)
df_good['file_des_1492R'] = df_good.apply(lambda x:os.path.join(folder_mnh,x['folder'],x['Moonid']+'-1492R.ab1'),axis=1)

df_good_wrong = df_good[~df_good.file_ori_1492R.apply(os.path.exists)]
df_good_wrong = df_good_wrong.copy()
df_good_wrong = df_good_wrong.iloc[:,0:3]
df_good_wrong.to_excel('/mnt/d/linux/W/moon/MoonZZD/nameNotes_fileNotExist.xlsx')

df_good_right = df_good[df_good.file_ori_1492R.apply(os.path.exists)]
df_good_right = df_good_right.copy()

import shutil
df_good_right.apply(lambda x:shutil.copy(x['file_ori_27F'],x['file_des_27F']), axis=1)
df_good_right.apply(lambda x:shutil.copy(x['file_ori_1492R'],x['file_des_1492R']), axis=1)

filename = '/mnt/d/linux/W/moon/MoonZZD/list.txt'
ls = open(filename).read().split()
folderSeqs = os.listdir('/mnt/d/linux/W/moon/MoonZZD/abiNew')
downloaded = [e.split('_')[0] for e in folderSeqs]


folder_abi = '/mnt/d/linux/W/moon/MoonZZD/abiNew/'
files = glob.glob(folder_abi + '**/*.ab1', recursive=True)
files_MNH = [e for e in files if 'MNH' in os.path.basename(e)]
for _f in files_MNH:
    _newfolder = os.path.dirname(_f).replace('abiNew','abi_MNH')
    if not os.path.exists(_newfolder):
        os.mkdir(_newfolder)
    shutil.copy(_f, _f.replace('abiNew','abi_MNH'))



## second part, HFY_ZGZ's data
import glob
folder_abi = '/mnt/d/linux/W/moon/MoonZZD/ab1HFY/'
files = glob.glob(folder_abi + '**/*.ab1', recursive=True)
print(len([e for e in files if '.ab1' in e]))
print(len([e for e in files if '-27F.ab1' in e]))
print(len([e for e in files if '-1492R.ab1' in e]))


import pandas as pd
df = pd.read_csv('/mnt/d/linux/W/moon/MoonZZD/nameNotesHFY.txt', sep='\t',dtype=str)

import os
folders = os.listdir(folder_abi)

folder_mnh = '/mnt/d/linux/W/moon/MoonZZD/ab1HFYnew/'
for _folder in folders:
    _folderpath = os.path.join(folder_mnh,_folder)
    if not os.path.exists(_folderpath):
        os.mkdir(_folderpath)

df_withna = df[df.Nameid.isna() | df.folder.isna()]
df_good = df.dropna()
df_good = df_good.copy()

df_good['file_ori_27F'] = df_good.apply(lambda x:os.path.join(folder_abi,x['folder'],x['Nameid']+'-27F.ab1'),axis=1)
df_good['file_ori_1492R'] = df_good.apply(lambda x:os.path.join(folder_abi,x['folder'],x['Nameid']+'-1492R.ab1'),axis=1)
df_good['file_des_27F'] = df_good.apply(lambda x:os.path.join(folder_mnh,x['folder'],x['Moonid']+'-27F.ab1'),axis=1)
df_good['file_des_1492R'] = df_good.apply(lambda x:os.path.join(folder_mnh,x['folder'],x['Moonid']+'-1492R.ab1'),axis=1)

df_good_wrong = df_good[~df_good.file_ori_1492R.apply(os.path.exists)]
df_good_wrong = df_good_wrong.copy()
df_good_wrong = df_good_wrong.iloc[:,0:3]
#df_good_wrong.to_excel('/mnt/d/linux/W/moon/MoonZZD/nameNotesHFY_fileNotExist.xlsx')

df_good_right = df_good[df_good.file_ori_1492R.apply(os.path.exists)]
df_good_right = df_good_right.copy()

import shutil
df_good_right.apply(lambda x:shutil.copy(x['file_ori_27F'],x['file_des_27F']), axis=1)
df_good_right.apply(lambda x:shutil.copy(x['file_ori_1492R'],x['file_des_1492R']), axis=1)


## clean data of Human Group
import glob
import os
from Bio import SeqIO
import io
from difflib import SequenceMatcher
import sys
sys.path.append('/mnt/c/Users/ATPs/Documents/GitHub/Moon/TaxaFinder/')
import abi
import itertools

folder = '/mnt/d/linux/W/moon/MoonZZD/abi_MNH/'
files = glob.glob(folder + '**/*.ab1', recursive=True)
## save files to fastq format


ls_fastq = []
for _f in files:
    _ab1 = SeqIO.read(_f,'abi')
    _ab1.id = os.path.basename(_f).split('.')[0]
    _fastq = _ab1.format("fastq")
    ls_fastq.append(_fastq)
open('/mnt/d/linux/W/moon/MoonZZD/20180814ZDD.fastq','w').write(''.join(ls_fastq))

## 20180815
ls_seqs = list(SeqIO.parse('/mnt/d/linux/W/moon/MoonZZD/20180814ZDD.fastq',format='fastq'))

ls_ids = [e.id.split('-')[0] for e in ls_seqs]
st_ids = set(ls_ids)
ls_ids = list(st_ids)
ls_idsTable = open('/mnt/d/linux/W/moon/MoonZZD/20180815MNH_idsWithSpecies.txt').read().split()
ls_idsTableMissing = [e for e in ls_idsTable if e not in st_ids]

## get merged sequences
dc_seqs = {e:[] for e in st_ids}
for _s in ls_seqs:
    dc_seqs[_s.id.split('-')[0]].append(_s)

from collections import Counter
print(Counter([len(e) for e in dc_seqs.values()]))

dc_seqsBest = {}
for _k, _v in dc_seqs.items():
    _LRcombinations = list(itertools.permutations(_v,2))
    _LRcombinationsGood = [[l,r] for l,r in _LRcombinations if '-27F' in l.id and '-1492R' in r.id]
    _seqs = [abi.mergeLRfastqStr(l,r) for l,r in _LRcombinationsGood]
    _seqsGood = [e for e in _seqs if e is not None]
    if len(_seqsGood) > 0:
        _seqBest = sorted(_seqsGood,key=len)[-1]
    dc_seqsBest[_k] = _seqBest

open('/mnt/d/linux/W/moon/MoonZZD/20180818MNH_ids_withAB1.txt','w').write('\n'.join(dc_seqs.keys()))
fout = open('/mnt/d/linux/W/moon/MoonZZD/20180818MNH_ids_withAB1.seqtab','w')
for _k,_v in dc_seqsBest.items():
    fout.write(_k+'\t'+_v+'\n')
fout.close()


#unzip all files
import zipfile
files_zip = glob.glob('/mnt/c/temp/**/*.zip', recursive=True)
for filename in files_zip:
    filezip = zipfile.ZipFile(filename)
    filezip.extractall(filename.replace('.zip',''))
#    break
files_seqAll = glob.glob('/mnt/c/temp/**/MNH*.seq', recursive=True)
print(len(files_seqAll))

import re
def cleanNNNsSeq(seq):
    '''
    seq is a DNA sequence with ATCGN(atcgn), return a longest seq with no n
    seq = 'nncgngtgcttacncatgcagtcgac', return 'catgcagtcgac'
    '''
    seqs = re.split(re.compile('n+',re.I), seq)
    seqs_sort = sorted(seqs, key=len)
    return seqs_sort[-1]

def getFastaFromFile(filename):
    '''
    given a filename, return a str of fasta format.
    some files contain extra info
    filename looks like '/mnt/d/linux/W/MoonNt/KC/180223 SEQ/MN210615.seq'
    return a string of fasta format
    '''
    fasta_str = ">" + re.split("\\\\|/",filename)[-1].split(".")[0] + "\n"
    txt = open(filename).read()
    if len(txt) == 0:
        print(filename,'empty file')
        return ''
    txt_ls = txt.split("\n")
    if txt_ls[0][0] != ">":
        seq = [e for e in txt_ls if re.match("[ATCGN]+ *\n?$",e,re.I)]
        if len(seq) == 0:#may contain RYSWKM, ambiguous bases
            seq = [e for e in txt_ls if re.match("[ATCGNRYSWKM]+ *\n?$",e,re.I)]
            if len(seq) == 0:
                print(filename, "No nucleotide sequences identified")
                return ""
            print(filename, "contains RYSWKM, ambiguous bases")
        fasta_seq = ''.join(seq)
        return fasta_str + cleanNNNsSeq(fasta_seq).upper() + "\n"
    seqs = list(SeqIO.parse(open(filename),format = "fasta"))
    if len(seqs) == 1:#if there is only one seq in fasta, use file name as sequence name
        fasta_seq = str(SeqIO.read(open(filename),format = "fasta").seq)
        return fasta_str + cleanNNNsSeq(fasta_seq).upper() + "\n"
    else:
        return "".join([">"+e.id+"\n"+cleanNNNsSeq(str(e.seq)).upper() +"\n" for e in seqs])

ls_seqAll =[getFastaFromFile(e) for e in files_seqAll]
ls_AllSeqIO = list(SeqIO.parse(io.StringIO(''.join(ls_seqAll)),'fasta'))
ls_mergedAllSeqIO = [e for e in ls_AllSeqIO if '27F' not in e.id and '1492R' not in e.id]
print(Counter([len(e.id) for e in ls_mergedAllSeqIO]))
#Counter({8: 4232, 9: 324, 10: 12, 16: 106, 18: 106})
dc_mergedAllSeqIO ={re.split('-|_',e.id)[0]:[] for e in ls_mergedAllSeqIO}


## 20180822 collect info for HumanMicrobiotaBank
### species identified in humans
