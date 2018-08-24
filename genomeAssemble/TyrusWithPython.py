#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 12:53:38 2018

@author: x
"""

# set up some constants
## program fastqc
FASTQC = '/home/p/FastQC/fastqc'
## program, Trim reads
TRIMMOMATIC = '/home/p/Trimmomatic-0.38/trimmomatic-0.38.jar'
## program, reads quality report and Trim
FASTP = '/home/p/fastp/fastp'
## program, Identify the best kmer for genome assembling
KMERGENIE = '/home/p/KmerGenie/kmergenie-1.7048/kmergenie'
## program, assemble the genome
ABYSS = '/home/p/abyss/bin/bin/abyss-pe'
## program, assemble the genome
VELVETH = '/home/p/anaconda/anaconda3_5.2.0/bin/velveth'
VELVETG = '/home/p/anaconda/anaconda3_5.2.0/bin/velvetg'
## program, assemble the genome
SPADES = '/home/p/anaconda/anaconda3_5.2.0/bin/spades.py'
## program, assemble the genome
FQ2FA = '/home/p/anaconda/anaconda3_5.2.0/bin/fq2fa'
IDBA_UD = '/home/p/anaconda/anaconda3_5.2.0/bin/idba_ud'
## program, assemble the genome
SOAPDENOVO = '/home/p/anaconda/anaconda3_5.2.0/bin/SOAPdenovo-127mer'
## program, merge assembled results
MERGE = '/home/p/CISA/CISA1.3/Merge.py'
CISA = '/home/p/CISA/CISA1.3/CISA.py'
## program, annotate merged genome
PROKKA = '/home/p/prokka/prokka-master/bin/prokka'



# set up some default values
THREADS = 24
OUTFOLDER = ''



import argparse
parser = argparse.ArgumentParser(description='Realization of Tychus with Python scripts. Add SOAPdenovo, improve the genus database for prokka with data from Mypro https://sourceforge.net/projects/sb2nhri/files/MyPro/ ')
