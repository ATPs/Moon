#!/usr/bin/env /mnt/d/linux/P/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 13:35:56 2018

@author: x
"""
folder = "."

import glob
import re
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import pickle


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

def getSeqFiles(folder):
    '''
    given a foler name, return a list of all sequence files.
    seq files are files end with .seq
    '''
    files = glob.glob(folder + '**/*.seq', recursive=True)
    files += glob.glob(folder + '**/*.SEQ', recursive=True)
    return files

def cleanNNNsSeq(seq):
    '''
    seq is a DNA sequence with ATCGN(atcgn), return a longest seq with no n
    seq = 'nncgngtgcttacncatgcagtcgac', return 'catgcagtcgac'
    '''
    seq = seq.upper()# change to capital letters
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

def taxaFinder(folder = folder, db = "fungi", taxa_order = "0"):
    '''
    given a path of a folder, output a file with the alignment and taxa info for the sequences
    db is the database to use. the keywords
    '''
    taxa = dcTaxa[taxa_order]
    # add the slash to make the path easier to use
    if folder[-1] != "/":
        folder = folder + "/"
    
    # get paths of all sequence files
    files = getSeqFiles(folder)
    
    # get all fasta sequences
    fasta_seqs = [getFastaFromFile(e) for e in files]
    file_prefix = folder + folder.split("/")[-2]+"_"
    # save the fasta_seqs to a file in folder
    open(file_prefix + "AllSeqs.fasta","w").write("".join(fasta_seqs))
    
    # run local blast against SILVA
    command_line = "/home/p/blast/bin/blastn -db {database} -query {folder}AllSeqs.fasta  -out {folder}AllSeqs.fasta.{keyword}.tab -num_threads 24 -max_target_seqs 100 -word_size 28 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp staxids'".format(folder = file_prefix, database = dcDB[db], keyword=db)
    subprocess.run(command_line, shell=True)
    
    # read in blast result for SILVA
    df_DB = pd.read_csv("{folder}AllSeqs.fasta.{keyword}.tab".format(folder = file_prefix, keyword=db), sep = "\t",header = None, dtype = {0:str,1:str,2:float,3:int,4:int,5:int,6:int,7:int,8:int,9:int,10:float,11:float,12:int,13:str})
    df_DB.columns = ["query", "subject","identity","matchLength","missMatch","gap","qstart","qend","sstart", "send", "evalue", "bitScore", "qcover", "taxID"]
    ## add column 14, length of identical bases
    df_DB["identical"] = round(df_DB["identity"] * df_DB["matchLength"]/100)
    
    
    df_DB["subject"]=df_DB["subject"].apply(lambda x:x.split("|")[1] if "|" in x else x)
    
    #for multiple match of a query/subject, only keept the first one
    df_DB = df_DB.drop_duplicates(subset=["query","subject"],keep='first') 
    
    # read in SILVA16sDNA.goodLineage as dataframe
    df_goodLineage = pickle.load(open('/mnt/d/linux/W/NCBI/taxonamy20180810/speciesLineageName.memory','rb'))
    
    df_DB = df_DB.join(df_goodLineage, on="taxID",how="left")
    df_DB.loc[:,['query', 'subject', 'matchLength', 'identical', 'qcover', 'identity', 'taxID','name','species_id'] + taxa].to_excel(file_prefix + "All{keyword}Result.xlsx".format(keyword=db), index = False)
    
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
    df_DB_keep1=df_DB_keep1.sort_index()
    
    dc_seq = SeqIO.to_dict(SeqIO.parse(open(file_prefix + "AllSeqs.fasta"),"fasta"))
    df_DB_keep1['seq'] = df_DB_keep1['query'].apply(lambda x:str(dc_seq[x].seq))
    df_DB_keep1['length'] = df_DB_keep1['seq'].apply(lambda x:len(x))
    
    df_DB_keep1.loc[:,['query', 'subject', 'seq', 'length', 'matchLength', 'identical', 'qcover','identity', 'taxID','name','species_id'] +  taxa].to_excel(file_prefix + "Best{keyword}Result.xlsx".format(keyword=db), index = False)
    
    queries = set(df_DB_keep1["query"])
    fout = open(file_prefix + "sequencesNotIdentified{keyword}.fasta".format(keyword=db),"w")
    for seq in dc_seq.values():
        if seq.id not in queries:
            fout.write(">"+seq.id+"\n"+str(seq.seq)+"\n")
    fout.close()

import argparse
parser = argparse.ArgumentParser(description='the folder path of sequencing results')
parser.add_argument('folder',help = 'input folder path of sequencing results')
parser.add_argument("--db", "-d", help = "database to use, default SILVA. Currently support SILVA, Ezbio", default = "SILVA")
parser.add_argument("--taxa_order", "-t", help = "taxa order. default 0, which is ['species', 'genus', 'family', 'order', 'class', 'phylum']", default = "0")
f = parser.parse_args()
folder = f.folder
if f.db in dcDB:
    taxaFinder(folder,db=f.db, taxa_order=f.taxa_order)
else:
    print("undefined database")