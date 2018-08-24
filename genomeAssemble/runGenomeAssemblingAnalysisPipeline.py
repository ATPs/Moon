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


folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene5species/raw/NE2_DMS26159-V/'
outprefix = 'Nitrosomonas_eutropha_18'
file_fq1 = 'NE2_1.fq'
file_fq2 = 'NE2_2.fq'
kmers = [59, 71, 83, 95, 107, 119, 131]
GenomeAssembling20180523.geneAssembling(folder, file_fq1, file_fq2, outprefix, kmers)


folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene5species/raw/NE4_DMS26160-V/'
outprefix = 'Nitrosomonas_eutropha_HJ-5'
file_fq1 = 'NE4_1.fq'
file_fq2 = 'NE4_2.fq'
kmers = [59, 71, 83, 95, 107, 119, 131]
GenomeAssembling20180523.geneAssembling(folder, file_fq1, file_fq2, outprefix, kmers)


folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene5species/raw/PP211_DMS26161-V/'
outprefix = 'MN213211'
file_fq1 = 'MN213211_1.fq'
file_fq2 = 'MN213211_2.fq'
kmers = [59, 71, 83, 95, 107, 119, 131]
GenomeAssembling20180523.geneAssembling(folder, file_fq1, file_fq2, outprefix, kmers)


folder = '/mnt/d/linux/W/MoonNGS/raw/20180515Novogene5species/raw/PPE15_DMS26162-V/'
outprefix = 'MN212471'
file_fq1 = 'MN212471_1.fq'
file_fq2 = 'MN212471_2.fq'
kmers = [59, 71, 83, 95, 107, 119, 131]
GenomeAssembling20180523.geneAssembling(folder, file_fq1, file_fq2, outprefix, kmers)


