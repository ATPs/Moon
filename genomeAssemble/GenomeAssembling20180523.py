#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 10:11:59 2018

@author: x
"""
# import modules
import os
import subprocess
from Bio import SeqIO
import glob
import itertools
import datetime
import shutil

# set up constants
#folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS1/'
#outprefix = 'Nitrosomonas_eutropha_JS1'
#file_fq1 = 'JS.1.L350_DMS22631-V_1.fq'
#file_fq2 = 'JS.1.L350_DMS22631-V_2.fq'
#kmers = [59, 71, 83, 95, 107, 119]


TXT_config = '''
#maximal read length
max_rd_len=1000
[LIB]
#average insert size
avg_ins=350
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
'''




TRIMMOMATIC = '/home/p/Trimmomatic-0.38/trimmomatic-0.38.jar'
SOAPDENOVO = '/home/p/SOAPdenovo2-master/SOAPdenovo-127mer'
GAPCLOSER = '/home/p/SOAPdenovo2-master/GapCloser'
ILLUMINACLIP = '/home/p/Trimmomatic-0.38/adapters/TruSeq3-PE.fa'
KMERFREQ = '/home/p/SOAPec_bin_v2.03/bin/KmerFreq_AR'
CORRECTOR = '/home/p/SOAPec_bin_v2.03/bin/Corrector_AR'
FASTQC = '/home/p/FastQC/fastqc'
file_config = 'soapdenovo_config.txt'
file_correct = 'soapdenovo_files.txt'
file_log = 'soapdenovo_logs.txt'

def TIME():
    a = datetime.datetime.now()
    t="%0d/%02d/%02d %02d:%02d:%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
    return t

def CountLenScaf(filename):
    '''
    counter number of scaffolds longer than 500bp
    '''
    n = 0
    for s in SeqIO.parse(open(filename),'fasta'):
        m = len(s.seq)
        if m > 500:
            n += m
    return n

def countNumScaf(filename):
    '''
    counter number of scaffolds longer than 500bp
    '''
    n = 0
    for s in SeqIO.parse(open(filename),'fasta'):
        if len(s.seq) > 500:
            n += 1
    return n

def countNumLenNsScaf(filename):
    '''
    count number of scaffolds longer than 500bp
    and the total length of scaffolds longer than 500bp
    and the number of Ns
    '''
    n = 0
    l = 0
    Ns = 0
    for s in SeqIO.parse(open(filename),'fasta'):
        m = len(s.seq)
        if m > 500:
            l += m
            n += 1
            Ns += str(s.seq).count('N')
    return n,l,Ns

def geneAssembling(folder, file_fq1, file_fq2, outprefix='', kmers=[35, 47, 59, 71, 83, 95, 107, 119]):
    # change to working directory
    os.chdir(folder)
    if not os.path.exists('soapdenovo'):
        os.mkdir('soapdenovo')
    os.chdir('soapdenovo')
    log_info = open(file_log,'a')
    log_info.write(TIME()+'  start run assembly\n')
    log_info.write('folder: '+folder+'\n')
    log_info.write('reads file 1: '+file_fq1+'\n')
    log_info.write('reads file 2: '+file_fq2+'\n')
    log_info.write('kmers: '+ ','.join(str(e) for e in kmers)+'\n\n')
    
    
    # check quality
    command_fastqc = '{FASTQC} ../{file_fq1} ../{file_fq2} -t 24'.format(FASTQC=FASTQC, file_fq1=file_fq1,file_fq2=file_fq2)
    subprocess.run(command_fastqc, shell=True)
    log_info.write('check fastq quality with fastqc\n\n')
    shutil.move('../{file_fq1}_fastqc.html'.format(file_fq1=file_fq1.rsplit('.',1)[0]),'./{file_fq1}_fastqc.html'.format(file_fq1=file_fq1.rsplit('.',1)[0]))
    shutil.move('../{file_fq2}_fastqc.html'.format(file_fq2=file_fq2.rsplit('.',1)[0]),'./{file_fq2}_fastqc.html'.format(file_fq2=file_fq2.rsplit('.',1)[0]))
    os.remove('../{file_fq1}_fastqc.zip'.format(file_fq1=file_fq1.rsplit('.',1)[0]))
    os.remove('../{file_fq2}_fastqc.zip'.format(file_fq2=file_fq2.rsplit('.',1)[0]))

    
    # Trimmomatic to trim reads
    command_trim = 'java -jar {TRIMMOMATIC}  PE -phred33 -threads 24 ../{file_fq1} ../{file_fq2}  {file_fq1}.pair  {file_fq1}.single {file_fq2}.pair {file_fq2}.single ILLUMINACLIP:{ILLUMINACLIP}:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:0 MINLEN:30'.format(TRIMMOMATIC=TRIMMOMATIC,  file_fq1=file_fq1, file_fq2=file_fq2, ILLUMINACLIP=ILLUMINACLIP)
    subprocess.run(command_trim, shell=True)
    log_info.write(TIME()+ '  Finished trim fastq reads with trimmoatic\n\n')
    
    # SOAPDENOVO

    with open(file_correct,'w') as f:
        f.write('{file_fq1}.pair\n'.format(file_fq1=file_fq1))
        f.write('{file_fq2}.pair\n'.format(file_fq2=file_fq2))
        f.write('{file_fq1}.single\n'.format(file_fq1=file_fq1))
        f.write('{file_fq2}.single\n'.format(file_fq2=file_fq2))
    log_info.write(TIME()+ '  Finished generate config file for soapdenovo\n\n')
    
    ## generate kmer counts to correct reads
    command_kmer = '{KMERFREQ}  -k 17 -t 24 -q 33 {file_correct}'.format(KMERFREQ=KMERFREQ, file_correct=file_correct)
    subprocess.run(command_kmer,shell=True)
    log_info.write(TIME()+ '  Finished counting kmers, kmer=17\n\n')
    
    ## Correct reads
    command_correct = '{CORRECTOR}  -k 17 -t 24 -Q 33 -o 3 output.freq.cz output.freq.cz.len {file_correct} '.format(CORRECTOR=CORRECTOR, file_correct=file_correct)
    subprocess.run(command_correct,shell=True)
    log_info.write(TIME()+ '  Finished correct reads, kmer=17\n\n')
    
    ## genverate config file
    ### don't know why some fastq sequences are in wrong format. Correct errors
    files_reads_single = []
    files_reads_single.append('{file_fq1}.single.cor.pair_1.fq'.format(file_fq1=file_fq1))
    files_reads_single.append('{file_fq2}.single.cor.pair_2.fq'.format(file_fq2=file_fq2))
    files_reads_single.append('{file_fq1}.pair.cor.single.fq'.format(file_fq1=file_fq1))
    files_reads_single.append('{file_fq1}.single.cor.single.fq'.format(file_fq1=file_fq1))
    fout = open('clean_single_reads.fq','w')
    for _f in files_reads_single:
        if os.path.isfile(_f):
            _ls = open(_f).readlines()
            for _l in range(0,len(_ls),4):
                if _ls[_l][0] == '@':
                    fout.write(_ls[_l] + _ls[_l+1] + _ls[_l+2] + _ls[_l+3])
    fout.close()
    
    with open(file_config,'w') as f:
        f.write(TXT_config)
        f.write('q1={file_fq1}.pair.cor.pair_1.fq\n'.format(file_fq1=file_fq1))
        f.write('q2={file_fq2}.pair.cor.pair_2.fq\n'.format(file_fq2=file_fq2))
#        f.write('q={file_fq1}.single.cor.pair_1.fq\n'.format(file_fq1=file_fq1))
#        f.write('q={file_fq2}.single.cor.pair_2.fq\n'.format(file_fq2=file_fq2))
#        f.write('q={file_fq1}.pair.cor.single.fq\n'.format(file_fq1=file_fq1))
#        f.write('q={file_fq1}.single.cor.single.fq\n'.format(file_fq1=file_fq1))
        f.write('q=clean_single_reads.fq\n')
    
    ## Generate commond lines
    command_soaps = []
    for kmer in kmers:
        command_soap = '{SOAPDENOVO} all -s {file_config} -K {kmer} -R -o {outprefix}.kmer{kmer} -d 1 -p 24'.format(SOAPDENOVO=SOAPDENOVO, folder=folder, file_config=file_config, kmer=kmer, outprefix=outprefix)
        command_soaps.append(command_soaps)
        subprocess.run(command_soap, shell=True)
        print(command_soap)
        log_info.write(TIME()+ '  Finished soapdenovo for kmer {kmer}\n'.format(kmer=kmer))
    
    ## select best kmer
    scaf_files = {kmer:'{outprefix}.kmer{kmer}.scafSeq'.format(kmer=kmer, outprefix=outprefix) for kmer in kmers}
    scaf_lens = {kmer:CountLenScaf(scaf_files[kmer]) for kmer in kmers}
    scaf_nums = {kmer:countNumScaf(scaf_files[kmer]) for kmer in kmers}
    mean_scaf_len = sum(scaf_lens.values()) / len(kmers)
    bestkmers = [kmer for kmer in kmers if scaf_lens[kmer] > mean_scaf_len]
    bestkmer = min(bestkmers,key=lambda x:scaf_nums[x])
    log_info.write('\n')
    log_info.write('summary for different kmer settings:\n')
    log_info.write('kmer\tscaf_num\tscaf_len\n')
    for kmer in kmers:
        log_info.write('{kmer}\t{sn}\t{sl}\n'.format(kmer=kmer, sn = scaf_nums[kmer], sl=scaf_lens[kmer]))
    log_info.write('bestkmer is '+str(bestkmer)+'\n\n')
    ## remove SOAPdenovo generated files
    files_rm = glob.glob(outprefix+".kmer*")
    for _f in files_rm:
        os.remove(_f)
    
    log_info.write('adjust -d -u -R -F to find best parameters\n')
    ## adjust -d -u -R -F to find best parameters
    lis_options = [' '.join(e) for e in itertools.product( ['-u', ''], ['-R', ''], ['-F', ''])]
    for _n in range(len(lis_options)):
        command_soap = '{SOAPDENOVO} all -s {file_config} -K {bestkmer} -o {outprefix}.bestkmer{_n} {option} -p 24'.format(SOAPDENOVO=SOAPDENOVO, folder=folder, file_config=file_config, bestkmer=bestkmer, outprefix=outprefix, _n=_n, option = lis_options[_n])
        subprocess.run(command_soap, shell=True)
        log_info.write(TIME()+ '  Finished soapdenovo for setting: {option}\n'.format(option=lis_options[_n]))
    
    for _n in range(len(lis_options)):
        command_gap = '{GAPCLOSER} -a {outprefix}.bestkmer{_n}.scafSeq -b {file_config} -o {outprefix}.bestkmer{_n}.scafSeqFilled -l 155 -t 24'.format(GAPCLOSER=GAPCLOSER, folder=folder, file_config=file_config, outprefix=outprefix, _n=_n)
        subprocess.run(command_gap, shell=True)
        log_info.write(TIME()+ '  Finished gapcloser for setting: {option}\n'.format(option = lis_options[_n]))
    
    scaf_files = ['{outprefix}.bestkmer{_n}.scafSeqFilled'.format(outprefix=outprefix,_n=_n) for _n in range(len(lis_options))]
    scaf_info = {_n: countNumLenNsScaf(scaf_files[_n]) for _n in range(len(lis_options))}
    
    min_scafnum = min(v[0] for v in scaf_info.values())
    min_scafnum_info = [(k,v) for k,v in scaf_info.items() if v[0]==min_scafnum]
    best_setting = max(min_scafnum_info, key=lambda x:x[1][1])
    best_option = lis_options[best_setting[0]]
    best_scaf_info = 'number of scaffold: {0}, length: {1}'.format(best_setting[1][0],best_setting[1][1])
    log_info.write('\nSummary for different soapdenovo options:\n')
    log_info.write('option:scaf_num\tscanf_len\tNs\n')
    for _n in range(len(lis_options)):
        log_info.write('"{option}":\t{sn}\t{sl}\t{Ns}\n'.format(option=lis_options[_n], sn = scaf_info[_n][0], sl = scaf_info[_n][1], Ns = scaf_info[_n][2]))
    log_info.write('\nBest option is: '+best_option)
    log_info.write('\nBest assembly has '+best_scaf_info+'\n\n')
    ## rename best assembly and remove files
    best_scaf_file = scaf_files[best_setting[0]]
    best_scaf_file_outname = outprefix+'.soapdenovo.scaffold.fasta'
    lis = list(SeqIO.parse(open(best_scaf_file),'fasta'))
    fout = open(best_scaf_file_outname,'w')
    for _s in lis:
        if len(_s.seq)>500:
            fout.write('>'+_s.description+'\n'+str(_s.seq)+'\n')
    fout.close()
    os.rename(best_scaf_file,outprefix+'.soapdenovo.scaffold.fasta')
    os.rename(outprefix+'.bestkmer'+str(best_setting[0])+'.scafStatistics', outprefix+'.soapdenovo.scaffold.info')
    files_rm = glob.glob(outprefix+".bestkmer*")
    for _f in files_rm:
        os.remove(_f)
    
    log_info.write(TIME()+" Finished")
    log_info.close()
    
    ## remove files
