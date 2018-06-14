#!/usr/bin/env /mnt/d/linux/P/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:25:40 2018

@author: x
"""

#testfile = '/mnt/d/linux/W/moon/Moon_CWJ/18-4/fungi/MN112848-ITS5.ab1'

from Bio import SeqIO
import re
import glob
import os

def cleanNNNsSeq(seq):
    '''
    seq is a DNA sequence with ATCGN(atcgn), return a longest seq with no n
    seq = 'nncgngtgcttacncatgcagtcgac', return 'catgcagtcgac'
    '''
    seqs = re.split(re.compile('n+',re.I), seq)
    seqs_sort = sorted(seqs, key=len)
    return seqs_sort[-1]

def readABI(filename):
    '''
    given a filename, read the file to SeqIO object
    '''
    f = open(filename,mode="rb")
    return SeqIO.read(f,format='abi')

def cleanABI(filename, qmin = 20, window = 10, minlen = 600, outputfile = False):
    '''
    given a filename of ABI format ".ab1", return a good quality sequence in fasta format.
    first read in the abi file, convert it to fastq format.
    then get the longest sequence with quality score >20. (average score in window)
    the score of Biopython is Phred+33
    if outputfile = True, output a fasta file.
    '''
    try:
        f = open(filename,mode="rb")
        s = SeqIO.read(f,format='abi')
        sq = s.format("fastq")
    except:
        print("file cannot be read", filename)
        return ""
    qmin += 33
    sqs = sq.split("\n")
    seq = sqs[1]
    quality = sqs[3]
    f.close()
    
    qualInt = list(map(ord,quality))
    maxlen = window -1
    start = -1
    end = -1
    selected = []
    for i in range(len(seq) - window + 1):
        q_mean = sum(qualInt[i:i+window]) / window
        if q_mean >= qmin:
            if "N" in seq[i:i+window]:# only works with N in the beginning of seq
#                print(i,"qmin:",qmin,"has N")
                continue
            if maxlen == window -1:
                start = i
                maxlen += 1
            else:
                maxlen += 1
        else:
            if start == -1:
                continue
            end = i+window - 1
#            print("maxlen:",maxlen)
            if maxlen >= window:
                selected.append((start, end, maxlen))
                maxlen = window -1
    if end < start + window -1:# scan till the end.
        end = len(seq)
        selected.append((start,end,end-start))
    
    if len(selected) < 1:
        print("The qualities of all bases are bad.")
        return ""
    
    bestseqs = [seq[e[0]:e[1]] for e in selected]
    bestseqs = [cleanNNNsSeq(e) for e in bestseqs]
    bestseq = sorted(bestseqs,key=len)[-1]
    
    if len(bestseq) <minlen:
        print("best sequence is shorter than",minlen,"for file",filename)
        return ""
    bestseq = ">"+sqs[0][1:]+"\n"+bestseq+"\n"
    if outputfile:
        open(filename.replace(".ab1",".seq"),"w").write(bestseq)
    return bestseq

def convertFolder2seq(folder,qmin = 15, window = 10, minlen = 300):
    '''
    convert .ab1 files in a folder to fasta formats
    '''
    files = glob.glob(os.path.join(folder,"**","*.ab1"), recursive=True)
    for filename in files:
        cleanABI(filename, qmin = qmin, window = window, minlen = minlen)