#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:25:35 2018

@author: x
"""
import os
import sys
dirname = os.path.dirname(__file__)
sys.path.append(dirname)

import GenomeAssembling20180523

#folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/JS4/'
#outprefix = 'Nitrosomonas_eutropha_JS4'
#file_fq1 = 'JS4V1.fq'
#file_fq2 = 'JS4V2.fq'
#kmers = [95,107, 119,131]

folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene3species/raw/ZD18/'
outprefix = 'Bacteroides_fragilis_ZD18'
file_fq1 = 'ZD18V_1.fq'
file_fq2 = 'ZD18V_2.fq'
kmers = [95, 107, 119,131]

GenomeAssembling20180523.geneAssembling(folder, file_fq1, file_fq2, outprefix, kmers)