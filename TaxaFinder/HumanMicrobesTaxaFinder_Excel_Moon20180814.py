#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 10:14:35 2018

@author: x
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 10:10:47 2018

@author: x
"""

#!/usr/bin/env /mnt/d/linux/P/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 13:35:56 2018

@author: x
"""

import glob
import re
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import pickle
import os
import io
import datetime
import shutil
import sys
import time

dcDB = {"fungi":"/mnt/d/linux/W/NCBI/fungi/ITS7.2",
        "SILVA":"/mnt/d/linux/W/NCBI/SILVA/SILVA16sDNA",
        "Ezbio":"/mnt/d/linux/W/NCBI/Ezbio20180810/Ezbio20180810",
        #"nt_micro":"/mnt/d/linux/W/NCBI/nt20180810/nt_micro",
        # use SSD to speedup
        "nt_micro":"/home/NCBI/nt_micro/nt_micro",
        }

DEFAULT = ['species', 'genus', 'family', 'order', 'class', 'phylum','kingdom']
REVERSE = ['kingdom','phylum', 'class', 'order', 'family', 'genus', 'species']
REVERSE2 = ['species', 'kingdom','phylum', 'class', 'order', 'family', 'genus']
dcTaxa = {"0":DEFAULT,"1":REVERSE,"2":REVERSE2}


def cleanNNNsSeq(seq):
    '''
    seq is a DNA sequence with ATCGN(atcgn), return a longest seq with no n
    seq = 'nncgngtgcttacncatgcagtcgac', return 'catgcagtcgac'
    '''
    seq = seq.upper()
    seqs = re.split(re.compile('n+',re.I), seq)
    seqs_sort = sorted(seqs, key=len)
    return seqs_sort[-1]
    
def getFastaFromExcelFile(df):
    '''
    given a  Excel filename, return a str of fasta format. Only with lines labeled "Y" in column "Analysis"
    ouput sequences for taxonomy
    and all available sequences
    '''
    df_analysis = df[df['Analysis']=='Y']
    df_all = df[df['seq'].notna()]
    seqs = []
    for _r, _e in df_analysis.iterrows():
        _seq = cleanNNNsSeq(_e['seq'])
        seqs.append('>' + _e['all_IDs'] +'\n' + _seq +'\n')
    seq = ''.join(seqs)
    seqs_all = []
    for _r, _e in df_all.iterrows():
        _seq = cleanNNNsSeq(_e['seq'])
        seqs_all.append('>' + _e['all_IDs'] +'\n' + _seq +'\n')
    seq_all = ''.join(seqs_all)
    return seq, seq_all

def check_species_name(name):
    '''
    check if species name is good.
    '''
    e = str(name)
    if e == 'nan':
        return False
    if re.search("environment", e, re.I) or re.search("unidentified", e, re.I) or re.search("uncultured", e, re.I):
        return False
    return True

def taxaFinder(filename,folder,outexcelname = 'MNH_human_micro_all',dbs=['Ezbio','nt_micro']):
    '''
    given a path of a folder, output a file with the alignment and taxa info for the sequences
    db is the database to use. the keywords
    '''
    #taxa = dcTaxa["0"]
    # add the slash to make the path easier to use
    
    
    #folder = '/mnt/d/linux/M/www/Django/moonbio/media/HumanMicrobes'
    #os.chdir(folder)
    currentfolder = folder
    a = datetime.datetime.now()
    foldernameYearMonthDay = "%0d%02d%02d"%(a.year,a.month,a.day)
    if os.path.exists(os.path.join(folder,foldernameYearMonthDay)):
        shutil.rmtree(os.path.join(folder,foldernameYearMonthDay))
    os.makedirs(os.path.join(folder,foldernameYearMonthDay))
    #os.chdir(foldernameYearMonthDay)
    
    # read in df_MNH excel file
    df_MNH = pd.read_excel(filename, sheet_name="all", index_col=0, dtype=str)
    df_MNH = df_MNH.replace('nan',np.nan)
    df_MNH = df_MNH.replace('',np.nan)
    
    # get all fasta sequences
    currentfolder = os.path.join(folder,foldernameYearMonthDay)
    file_prefix = os.path.join(currentfolder,foldernameYearMonthDay)
    fasta_seqs, fasta_seqs_all = getFastaFromExcelFile(df_MNH)
    open(file_prefix+'query.fasta','w').write(fasta_seqs)
    open(file_prefix+'all.fasta','w').write(fasta_seqs_all)
    
    
    # run local blast against Ezbio
    db = dbs[0]
    command_line = "/home/p/blast/bin/blastn -db {database} -query {folder}query.fasta  -out {folder}query.fasta.{keyword}.tab -num_threads 24 -max_target_seqs 100 -word_size 28 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp staxids'".format(folder = file_prefix, database = dcDB[db], keyword=db)
    print(command_line)
    #subprocess.run(command_line, shell=True)
    os.system(command_line)
    time.sleep(5)
    df_DB1 = pd.read_csv("{folder}query.fasta.{keyword}.tab".format(folder = file_prefix, keyword=db), sep = "\t",header = None, dtype = {0:str,1:str,2:float,3:int,4:int,5:int,6:int,7:int,8:int,9:int,10:float,11:float,12:int,13:str})
    df_DB1.columns = ["query", "subject","identity","matchLength","missMatch","gap","qstart","qend","sstart", "send", "evalue", "bitScore", "qcover", "taxID"]
    print('blast Ezbio shape', df_DB1.shape)
    #for multiple match of a query/subject, only keept the first one
    df_DB1 = df_DB1.drop_duplicates(subset=["query","subject"],keep='first')
    df_DB1["identical"] = round(df_DB1["identity"] * df_DB1["matchLength"]/100)
    df_DB1["subject"]=df_DB1["subject"].apply(lambda x:x.split("|")[1] if "|" in x else x)
    print('blast Ezbio shape, remove duplicated subject for same query', df_DB1.shape)
    
    db = dbs[1]
    command_line = "/home/p/blast/bin/blastn -db {database} -query {folder}query.fasta  -out {folder}query.fasta.{keyword}.tab -num_threads 24 -max_target_seqs 100 -word_size 28 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp staxids'".format(folder = file_prefix, database = dcDB[db], keyword=db)
    print(command_line)
    #subprocess.run(command_line, shell=True)
    os.system(command_line)
    time.sleep(5)
    
    df_DB2 = pd.read_csv("{folder}query.fasta.{keyword}.tab".format(folder = file_prefix, keyword=db), sep = "\t",header = None, dtype = {0:str,1:str,2:float,3:int,4:int,5:int,6:int,7:int,8:int,9:int,10:float,11:float,12:int,13:str})
    df_DB2.columns = ["query", "subject","identity","matchLength","missMatch","gap","qstart","qend","sstart", "send", "evalue", "bitScore", "qcover", "taxID"]
    print('blast nt_micro shape', df_DB2.shape)
    df_DB2 = df_DB2.drop_duplicates(subset=["query","subject"],keep='first')
    print('blast nt_micro shape, remove duplicated subject for same query', df_DB2.shape)
    ## add column 14, length of identical bases
    df_DB2["identical"] = round(df_DB2["identity"] * df_DB2["matchLength"]/100)
    df_DB2["subject"]=df_DB2["subject"].apply(lambda x:x.split("|")[1] if "|" in x else x)
    
    
    
    # read in lineage
    df_goodLineage_hm = pd.read_csv('/mnt/d/linux/W/NCBI/MoonRef20180810/HumanRef2018.speciesLineageName',dtype=str,sep='\t')
    df_goodLineage_hm.index = df_goodLineage_hm['tax_id']
    df_goodLineage_hm = df_goodLineage_hm.drop(labels=['tax_id'],axis='columns')
    # read in ncbi taxa as dataframe
    df_goodLineage = pickle.load(open('/mnt/d/linux/W/NCBI/taxonamy20180810/speciesLineageName.memory','rb'))
    print("finished reading taxa lineage")
    
    # find human microbial species
    df_DB_hm1 = df_DB1.copy()
    df_DB_hm2 = df_DB2.copy()
    df_DB_hm1 = df_DB_hm1.join(df_goodLineage, on="taxID",how="left")
    df_DB_hm2 = df_DB_hm2.join(df_goodLineage_hm, on="taxID",how="left")
    df_DB_hm = pd.concat([df_DB_hm1,df_DB_hm2],axis='index',ignore_index=True)
    df_DB_hm = df_DB_hm.sort_values(by = ["query","bitScore","identity"], ascending=[True,False,False])
    df_DB_hm['goodLineage'] = df_DB_hm.apply(lambda x:x.isnull().sum() < 3 and check_species_name(x['species']) ,axis=1)
    df_DB_hm_1 = df_DB_hm[df_DB_hm['goodLineage']]
    df_DB_1 = df_DB_hm_1.copy()
    good_bitscore = [True]
    for _n in range(df_DB_1.shape[0]-1):
        if df_DB_1.iloc[_n,0] == df_DB_1.iloc[_n+1,0]:
            if df_DB_1.iloc[_n,2] == df_DB_1.iloc[_n+1,2] and df_DB_1.iloc[_n,11] == df_DB_1.iloc[_n+1,11]:
                good_bitscore.append(good_bitscore[_n])
            else:
                good_bitscore.append(False)
        else:
            good_bitscore.append(True)
    df_DB_1['goodScore']=good_bitscore
    df_DB_goodbit = df_DB_1[good_bitscore].copy()
    df_DB_goodbit['goodSpecies'] = df_DB_goodbit['species'].apply(lambda x:' sp.' not in x)
    df_DB_goodbit=df_DB_goodbit.sort_values(by=['query','goodSpecies'],ascending=False)
    df_DB_keep1 = df_DB_goodbit.drop_duplicates(subset = ["query"], keep = "first").copy()
    df_DB_hm_final = df_DB_keep1.sort_index()
    print('df_DB_hm_final shape', df_DB_hm_final.shape)
    
    
    ## find nt microbial species
    df_DB = pd.concat([df_DB1,df_DB2],axis='index',ignore_index=True)
    df_DB = df_DB.sort_values(by = ["query","bitScore","identity"], ascending=[True,False,False])
    
    df_DB = df_DB.join(df_goodLineage, on="taxID",how="left")
    df_DB['goodLineage'] = df_DB.apply(lambda x:x.isnull().sum() < 3 and check_species_name(x['species']) ,axis=1)
    df_DB_1 = df_DB[df_DB['goodLineage']].copy()
    good_bitscore = [True]
    for _n in range(df_DB_1.shape[0]-1):
        if df_DB_1.iloc[_n,0] == df_DB_1.iloc[_n+1,0]:
            if df_DB_1.iloc[_n,2] == df_DB_1.iloc[_n+1,2] and df_DB_1.iloc[_n,11] == df_DB_1.iloc[_n+1,11]:
                good_bitscore.append(good_bitscore[_n])
            else:
                good_bitscore.append(False)
        else:
            good_bitscore.append(True)
    df_DB_1['goodScore']=good_bitscore
    df_DB_goodbit = df_DB_1[good_bitscore].copy()
    df_DB_goodbit['goodSpecies'] = df_DB_goodbit['species'].apply(lambda x:' sp.' not in x)
    df_DB_goodbit=df_DB_goodbit.sort_values(by=['query','goodSpecies'],ascending=False)
    df_DB_keep1 = df_DB_goodbit.drop_duplicates(subset = ["query"], keep = "first").copy()
    df_DB_nt_final=df_DB_keep1.sort_index()
    print('df_DB_nt_final shape', df_DB_nt_final.shape)
    
    # merge the info in df_DB_hm_final and df_DB_nt_final to df_MNH
    mnh_index = df_MNH.index
    df_MNH.index = df_MNH['all_IDs']
    df_DB_hm_final.index = df_DB_hm_final['query']
    df_DB_nt_final.index = df_DB_nt_final['query']
    cols = ['subject', 'identity', 'matchLength', 'qcover', 'taxID', 'identical', 'species_id', 'name', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    for _r in df_DB_hm_final.index:
        df_MNH.loc[_r,'Analysis'] = np.nan
        for _c in cols:
            df_MNH.loc[_r,'HM_'+_c] = df_DB_hm_final.loc[_r,_c]
    for _r in df_DB_nt_final.index:
        for _c in cols:
            df_MNH.loc[_r,'nt_'+_c] = df_DB_nt_final.loc[_r,_c]
    print('finish update df_MNH with df_DB_hm_final and df_DB_nt_final')
    
    # blast to MNH DB
    sys.path.append('/mnt/d/linux/M/www/Django/TaxaFinder')
    import largeFastaFile
    ls_MNH = largeFastaFile.open_fasta_to_list(file_prefix+'all.fasta')
    dc_MNHLen = {e.id:len(e.seq) for e in ls_MNH}
    kmerlen = min(len(e) for e in ls_MNH) - 20
    ls_MNHunique = largeFastaFile.fasta_within_seq_big_withError(ls_MNH,error_rate = 0,kmerlen = kmerlen)
    files_toRemove = glob.glob(folder+'/moonRef/*.*')
    hm_ref_file = folder+'/moonRef/'+foldernameYearMonthDay+'HM_ref'
    largeFastaFile.saveFastaListToFile(ls_MNHunique,hm_ref_file)
    print('finish find non-redundant sequences, before',len(ls_MNH),'after',len(ls_MNHunique))
    
    
    os.system('/home/p/blast/bin/makeblastdb -in {hm_ref_file} -input_type fasta -dbtype nucl -out {hm_ref_file} -title "hm_ref" -parse_seqids'.format(hm_ref_file=hm_ref_file))
    command_line = "/home/p/blast/bin/blastn  -query {folder}all.fasta -db {hm_ref_file}  -out {folder}all.fasta.{keyword}.tab -num_threads 24 -max_target_seqs 1 -word_size 28 -outfmt '6' ".format(folder = file_prefix, hm_ref_file=hm_ref_file, keyword='hm_ref')
    subprocess.run(command_line, shell=True)
    df_DB_hm_ref = pd.read_csv('{folder}all.fasta.hm_ref.tab'.format(folder = file_prefix), sep='\t',header = None, dtype=str)
    dc_MNH2ref = dict(zip(df_DB_hm_ref[0], df_DB_hm_ref[1]))
    [os.remove(e) for e in files_toRemove if foldernameYearMonthDay not in e]
    for _k,_v in dc_MNH2ref.items():
        df_MNH.loc[_k,'MatchMNH'] = _v
        df_MNH.loc[_k,'MatchMNH_len'] = dc_MNHLen[_k]
    print('finish update match hm_ref_MNH')
    
    #save the result
    df_MNH.index = mnh_index
    df_MNH.to_excel(os.path.join(folder, foldernameYearMonthDay, foldernameYearMonthDay+outexcelname+'.xlsx'), sheet_name='all')
    ## move the output and clean the intermediate files
    os.rename(os.path.join(folder, foldernameYearMonthDay, foldernameYearMonthDay+outexcelname+'.xlsx'), os.path.join(folder, foldernameYearMonthDay+outexcelname+'.xlsx'))
    os.rename(file_prefix+'all.fasta',os.path.join(folder, 'moonAll',foldernameYearMonthDay+'all.fasta'))
    print("finish writing the final result")
    
    #os.chdir('..')
    shutil.rmtree(os.path.join(folder, foldernameYearMonthDay))
    return foldernameYearMonthDay

#filename = '/mnt/d/linux/M/www/Django/moonbio/media/HumanMicrobes/20180816MNH_human_micro_all_2.xlsx'
#taxaFinder(filename)