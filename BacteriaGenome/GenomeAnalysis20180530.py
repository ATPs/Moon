#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 09:41:41 2018

@author: x
"""
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import glob
import itertools
import datetime
import shutil

# constant values
QUAST = '/mnt/d/linux/P/quast/quast-4.6.3/quast.py'
GENEMARK = '/mnt/d/linux/P/quast/quast-4.6.3/quast_libs/genemark/linux_64/gmsn.pl'
BLASTP = '/mnt/d/linux/P/blast/bin/blastp'
VFDB_PR = '/mnt/d/linux/W/NCBI/VFDB/VFDB_blast201806/VFDB_Pr'
VFDB_INFO = '/mnt/d/linux/W/NCBI/VFDB/VFDB_info.txt'
VFDB_HEAD = '/mnt/d/linux/W/NCBI/VFDB/VFDB_blast201806/VFDB_Pr.head'
VFDB_VF = '/mnt/d/linux/W/NCBI/VFDB/VFDB_blast201806/VFDB_Pr.VF'
KO_info = '/mnt/d/linux/W/NCBI/KEGG/KO.txt'

# variables
REF_GENOME = '/mnt/d/linux/W/species/Nitrosomonas_eutropha/C91/GCF_000014765.1_ASM1476v1/GCF_000014765.1_ASM1476v1_genomic.fna'
REF_GENE = '/mnt/d/linux/W/species/Nitrosomonas_eutropha/C91/GCF_000014765.1_ASM1476v1/GCF_000014765.1_ASM1476v1_genomic.gff'
test_genome = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS1/soapdenovo/Nitrosomonas_eutropha_JS1.soapdenovo.scaffold.fasta'
outfolder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS1/soapdenovo/'

REF_GENOME = '/mnt/d/linux/W/species/Nitrosomonas_eutropha/C91/GCF_000014765.1_ASM1476v1/GCF_000014765.1_ASM1476v1_genomic.fna'
REF_GENE = '/mnt/d/linux/W/species/Nitrosomonas_eutropha/C91/GCF_000014765.1_ASM1476v1/GCF_000014765.1_ASM1476v1_genomic.gff'
test_genome = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS4/soapdenovo/Nitrosomonas_eutropha_JS4.soapdenovo.scaffold.fasta'
outfolder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS4/soapdenovo/'

REF_GENOME = '/mnt/d/linux/W/species/Bacteroides_fragilis/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.fna'
REF_GENE = '/mnt/d/linux/W/species/Bacteroides_fragilis/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.gff'
test_genome = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/ZD18/soapdenovo/Bacteroides_fragilis_ZD18.soapdenovo.scaffold.fasta'
outfolder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/ZD18/soapdenovo/'

def gene_mark(test_genome=test_genome,outfolder=outfolder):
    '''
    predict genes
    '''
    folder = os.path.join(outfolder,'genemark')
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    command_genemark = '{GENEMARK} --format GFF --clean --gm --fnn --faa {test_genome}'.format(GENEMARK=GENEMARK, test_genome=test_genome)
    subprocess.run(command_genemark,shell=True)
    file_aa = os.path.basename(test_genome)+'.faa'
    file_nn = os.path.basename(test_genome)+'.fnn'
    
    def rename(file):
        lis = open(file).readlines()
        fout = open(file,'w')
        for e in lis:
            if len(e)>0:
                if e[0]=='>':
                    e = e.replace('|',' ',1)
                fout.write(e)
        fout.close()
    
    rename(file_aa)
    rename(file_nn)
    

def VFDB_analysis(outfolder=outfolder, test_genome=test_genome, VFDB_score_min = 0.6):
    
    folder = os.path.join(outfolder,'VFDB')
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    file_pr = os.path.basename(test_genome) + '.faa'
    command = '{BLASTP} -query {outfolder}/genemark/{file_pr} -db {VFDB_PR} -outfmt 6 -num_threads 4 -out result_VFDB.tab -max_target_seqs 100 -task blastp-fast '.format(BLASTP=BLASTP, outfolder=outfolder, file_pr=file_pr, VFDB_PR=VFDB_PR)
    subprocess.run(command, shell=True)
    dc_ref = SeqIO.to_dict(SeqIO.parse(open(VFDB_PR),'fasta'))
    dc_ref_len = {k:len(v) for k,v in dc_ref.items()}
    dc_test = SeqIO.to_dict(SeqIO.parse(open(outfolder+'/genemark/'+file_pr),'fasta'))
    dc_test_len = {k:len(v) for k,v in dc_test.items()}
    df_align = pd.read_csv('result_VFDB.tab', sep='\t',header=None, names=["query", "subject","identity","matchLength","missMatch","gap","qstart","qend","sstart", "send", "evalue", "bitScore"])
    df_align_select = df_align[df_align.loc[:,'identity']>70]
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
    
    
    

def quast_analysis(REF_GENOME=REF_GENOME,REF_GENE=REF_GENE, test_genome=test_genome, outfolder=outfolder):
    '''
    
    '''
    
    command_quast = '{QUAST} -R {REF_GENOME} -G {REF_GENE}  {test_genome} -o {outfolder}/quast -t 4'.format(QUAST=QUAST, REF_GENOME=REF_GENOME, REF_GENE=REF_GENE, test_genome=test_genome, outfolder=outfolder)
    subprocess.run(command_quast, shell=True)

def kegg_generateComparison(file_test, file_ref):
    '''
    use https://www.kegg.jp/blastkoala/ to get the KOs of translated proteins
    extract those with KO and compare with reference KO list
    '''
    file_test = '/media/sf_linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS1/soapdenovo/KEGG/user_ko.txt'
    file_ref = '/media/sf_linux/W/species/Nitrosomonas_eutropha/kegg/ko_net.txt'
    
    file_test = '/media/sf_linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS4/soapdenovo/KEGG/user_ko.txt'
    file_ref = '/media/sf_linux/W/species/Nitrosomonas_eutropha/kegg/ko_net.txt'
    
    file_test = '/media/sf_linux/W/MoonNGS/raw/20180515Novogene3species/raw/ZD18/soapdenovo/KEGG/user_ko.txt'
    file_ref = '/media/sf_linux/W/species/Bacteroides_fragilis/GCF_000009925.1_ASM992v1/ko.txt'
    
    
    df = pd.read_csv(file_ref,sep='\t',header=None,names=['ref_id','ko'])
    df['ko'] = df['ko'].apply(lambda x:x.split(':')[-1])
    df['ref_id'] = df['ref_id'].apply(lambda x:x.split(':')[-1])
    df_t = pd.read_csv(file_test,sep='\t',header=None,names=['test_id','ko'])
    df_t = df_t.dropna()
    dfk = df.groupby('ko').agg(';'.join)
    df_tk = df_t.groupby('ko').agg(';'.join)
    df_f = dfk.join(df_tk,how='outer')
    dfko = pd.read_csv(KO_info,sep='\t',header=None,names=['ko','ko_info'], index_col=0)
    df_f = df_f.join(dfko,how='left')
    outfile = os.path.join(os.path.dirname(file_test),'ko_compare_with_ref.xlsx')
    df_f.to_excel(outfile)
    outfile2 = os.path.join(os.path.dirname(file_test),'ko_compare_with_ref.txt')
    fout = open(outfile2,'w')
    fout.write('#test_genome\n')
    for line in df_t.iterrows():
        fout.write('\t'.join(line[1]) +'\n')
    fout.write('#ref_genome\n')
    for line in df.iterrows():
        fout.write('\t'.join(line[1]) +'\n')
    fout.close()

def kegg_analysis():
    '''
    reconstruct pathway https://www.kegg.jp/kegg/tool/map_pathway.html
    download the webpage, and use the program below to download figures
    '''
    import re
    import requests
    import urllib
    folder = os.path.join(outfolder,'KEGG')
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    
    file_map = '/media/sf_linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS1/soapdenovo/KEGG/KEGG_map.html'
    file_map = '/media/sf_linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS4/soapdenovo/KEGG/Reconstruct PATHWAY.html'
    file_map = '/media/sf_linux/W/MoonNGS/raw/20180515Novogene3species/raw/ZD18/soapdenovo/KEGG/Reconstruct PATHWAY.html'
    
    response = open(file_map).read()
    pattern = re.compile('href="https://www.kegg.jp/kegg-bin/show_pathway\?\d*/map\d*\.coords\+reference"')
    keggs = re.findall(pattern,response)
    print(len(keggs))
    keggs = [e[6:-1] for e in keggs]
    image = re.compile('src="/tmp/mark_pathway\d*/map\d*.*?\.png"')
    metabolism_maps = '''00100 00270 00350 00471 00550 00630 00750 00970 01230 03010 04122 00010 00121 00280 00360 00472 00561 00633 00760 01040 01501 03018 00020 00130 00281 00362 00473 00562 00640 00770 01100 01502 03020 00030 00190 00290 00380 00480 00564 00650 00780 01110 01503 03030 00040 00220 00300 00400 00500 00590 00660 00785 01120 02010 03060 00051 00230 00310 00401 00520 00592 00670 00790 01130 02020 03070 00052 00240 00311 00410 00521 00620 00680 00860 01200 02024 03410 00053 00250 00330 00430 00523 00623 00710 00909 01210 02030 03420 00061 00260 00332 00450 00525 00625 00730 00910 01212 02040 03430 00071 00261 00340 00460 00540 00626 00740 00920 01220 02060 03440'''.split()
    keggs_keep = [e for e in keggs if e[-22:-17] in metabolism_maps]
    for _k in keggs_keep:
        _t = requests.get(_k).text
        imageurl = re.findall(image,_t)
        if len(imageurl)!=1:
            print(_k, imageurl)
        else:
            imageurl = imageurl[0]
            imageurl = 'https://www.kegg.jp/'+imageurl[5:-1]
            imagename = re.findall('map\d*',imageurl)[0]+'.png'
            urllib.request.urlretrieve(imageurl,imagename)
