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
from difflib import SequenceMatcher
import io

def cleanNNNsSeq(seq):
    '''
    seq is a DNA sequence with ATCGN(atcgn), return a longest seq with no n
    seq = 'nncgngtgcttacncatgcagtcgac', return 'catgcagtcgac'
    '''
    seqs = re.split(re.compile('n+',re.I), seq)
    seqs_sort = sorted(seqs, key=len)
    return seqs_sort[-1]

def readABI2fastq(filename):
    '''
    given a filename, read the file to SeqIO object
    '''
    if type(filename) is SeqIO.SeqRecord:
        return filename.format('fastq')
    try:
        f = open(filename,mode="rb")
        s = SeqIO.read(f,format='abi')
        sq = s.format("fastq")
        f.close()
        return sq
    except:
        print("file cannot be read", filename)
        return None

def fastqFilterWithWindow(sq, qmin = 20, window = 10):
    '''
    given a fastq sequence in txt format.
    '@id\nATCG\n+\nQQQQ\n'
    return a fastq seq with desired quality
    '''
    sqs = sq.split("\n")
    seqid = sqs[0][1:]
    seq = sqs[1]
    quality = sqs[3]
    qualInt = list(map(ord,quality))
    qualInt = [e -33 for e in qualInt]
    
    #find the longest part with no N
    seq_l = re.findall('[ATCG]+',seq)
    if len(seq_l) < 1:
        print(seqid, "empty!!!")
        return None
    seq_l = sorted(seq_l,key=len)[-1]
    seq_l_start = seq.find(seq_l)
    qualInt_l = qualInt[seq_l_start:seq_l_start+len(seq_l)]
    qual_l = quality[seq_l_start:seq_l_start+len(seq_l)]
    
    q_mins = [sum(qualInt_l[i:i+window]) / window for i in range(len(seq_l) - window + 1)]
    _n = 0
    good_regions = []
    while _n < len(q_mins):
        _q = q_mins[_n]
        if _q >qmin:
            _m = _n+1
            while _m < len(q_mins):
                _q = q_mins[_m]
                if _q > qmin:
                    _m += 1
                else:
                    good_regions.append([_n,_m])
                    _n = _m + 1
                    break
            else:
                good_regions.append([_n,_m])
                _n += 1
        else:
            _n += 1
    if len(good_regions) == 0:
        print("nogood")
        return None
    good_seq_regions = [[m,n+10-1] for m,n in good_regions]
    best_region = sorted(good_seq_regions,key=lambda x:x[1]-x[0])[-1]
    best_seq = seq_l[best_region[0]:best_region[1]]
    best_qual = qual_l[best_region[0]:best_region[1]]
    return '@' + seqid +'\n'+best_seq+'\n'+'+'+best_qual+'\n'
    
    
    
    
def cleanABI(filename, qmin = 20, window = 10, minlen = 200, outputfile = False):
    '''
    given a filename of ABI format ".ab1", return a good quality sequence in fasta format. filename can also be a SeqIO.SeqRecord
    first read in the abi file, convert it to fastq format.
    then get the longest sequence with quality score >20. (average score in window)
    the score of Biopython is Phred+33
    if outputfile = True, output a fasta file.
    '''
    sq = readABI2fastq(filename)
    if sq is None:
        return ""
    
    bestseq_fastq = fastqFilterWithWindow(sq=sq,qmin=qmin,window=window)
    if bestseq_fastq is None:
        print("best sequence is shorter than",minlen,"for file",filename)
        return ""
    sqs = bestseq_fastq.split('\n')
    bestseq = sqs[1]
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

def merge_seqStr(seql_seq, seqrc_seq, min_identical = 50):
    '''
    seql_seq and seqrc_seq are strings of sequence 
    merge them if their identical bases are longer than min_identical
    else, return the logest seql_seq or seqrc_seq
    '''
    seq_match = SequenceMatcher(None, seql_seq, seqrc_seq,autojunk=False).find_longest_match(0, len(seql_seq), 0, len(seqrc_seq))
    if seq_match.size >= min_identical:
        return seql_seq[:seq_match.a+seq_match.size] + seqrc_seq[seq_match.b+seq_match.size:]
    if len(seql_seq) > len(seqrc_seq):
        return seql_seq
    return seqrc_seq

def mergeLRfastqStr(seqlq, seqrq, min_identical = 50):
    '''
    seqlq and seqrq are seqIO in fastq format
    clean them first. Then merge.
    Return a string of sequence
    '''
    seqr = cleanABI(seqrq,minlen=1)
    seql = cleanABI(seqlq,minlen=1)
    if seqr == "":
        if seql == "":
            return ""
        return seql.split()[1]
    if seql == "":
        if seqr == "":
            return ""
        return seqr.split()[1]
    seql_seq = seql.split()[1]
    seqrc_seq = str(SeqIO.read(io.StringIO(seqr),'fasta').reverse_complement().seq)
    seq_merge = merge_seqStr(seql_seq=seql_seq, seqrc_seq=seqrc_seq, min_identical = min_identical)
    return seq_merge

def mergeLRfastaStr(seql, seqr, min_identical = 50):
    '''
    seqlq and seqrq are seqIO in fasta format
    Return a string of sequence
    '''
    if str(seqr.seq) == "":
        return str(seql.seq)
    if str(seql.seq) == "":
        return str(seqr.seq)
    seqrc_seq = str(seqr.reverse_complement().seq)
    seq_merge = merge_seqStr(seql_seq=str(seql.seq), seqrc_seq=seqrc_seq, min_identical = min_identical)
    return seq_merge

def mergerFastaSeqInFile(filename):
    '''
    given a filename of fasta sequences, return a str of merged sequences, and save to the original file. Also, give some discription about no good sequences.
    '''
    #filename = '/mnt/d/linux/M/2018081609305280-199244728_ab1.seq'
    
    lis = list(SeqIO.parse(open(filename),'fasta'))
    dc = {}
    for e in lis:
        _k = e.id.split('-')[0]
        if _k not in dc:
            dc[_k] = []
        dc[_k].append(e)
    
    fout = open(filename,'w')
    
    message = '\n'
    seq2 = ''
    
    for _k,_v in dc.items():
        if not re.match('^MNH?\d{5,10}$',_k):
            message = message + _k+': naming wrong.\n'
            continue
        if len(_v) == 1:
            if len(_v[0].seq) < 600:
                message = message + _k +': one end. Shorter than 600bp\n'
            else:
                _s = '>'+_k+'\n'+str(_v[0].seq)+'\n'
                fout.write(_s)
                seq2 += _s
        elif len(_v) == 2 and _v[0].id in [_k+'-1492R',_k+'-27F'] and _v[1].id in [_k+'-1492R',_k+'-27F']:
            _seq = mergeLRfastaStr(_v[0], _v[1])
            if len(_seq) < 600:
                message = message + _k +':two end. Shorter than 600bp\n'
            else:
                _s = '>'+_k+'\n'+_seq+'\n'
                fout.write(_s)
                seq2 += _s
        else:
            message = message + _k +': wrong name or more than 2 seqs under the ID\n'
    
    fout.close()
    print(message)
    return message, seq2
            

