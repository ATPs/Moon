filename = "/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018_from_XianzhiJiang.txt"

import pandas as pd
df_gut = pd.read_csv(filename,sep='\t')

file_name = "/mnt/d/linux/W/NCBI/taxonamy/names.dmp"
df_name = pd.read_csv(file_name,sep = "\t", header = None,dtype = str)
df_name = df_name.loc[:,[0,2,6]]
df_name.index = df_name[2]
df_gut2 = df_gut.join(df_name, on='Name',how='left')
df_gut.shape
df_gut2.shape
from collections import Counter
df_gut2 = df_gut2.drop_duplicates(subset = "Name")

df_gut3 = df_gut2[df_gut2[0].notnull()]
df_gut4 = df_gut3.iloc[:,0:3]


df_gut5 = df_gut2[df_gut2[0].isnull()]
txt = open(file_name).read()
[e for e in df_gut5.Name if e in txt]
txt = txt.upper()
[e for e in df_gut5.Name if e.upper() in txt]




## find 16s sequences of 2690 species that is identified in human
### SILVA some taxID are subspecies, change to species level
file_node = "/mnt/d/linux/W/NCBI/taxonamy/nodes.dmp"
import pandas as pd
df_node = pd.read_csv(file_node, sep = "\t",header = None,dtype = str)
df_node = df_node.loc[:,[0,2,4]]
df_node.columns = ["From","To","rank"]
df_node.head()
set(df_node.loc[:,"rank"])
from collections import Counter
print(Counter(df_node["rank"]))
dc_node = dict(zip(list(df_node.loc[:,"From"]), list(df_node.loc[:,"To"])))
dc_rank = dict(zip(list(df_node.loc[:,"From"]), list(df_node.loc[:,"rank"])))

### read in tax ids
file = "/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018.txt"
import pandas as pd
df_species = pd.read_csv(file, sep='\t',dtype=str)
taxids = list(df_species.taxid)
taxids_s = set(taxids)
## read in SILVA 16s
file = "/mnt/d/linux/W/NCBI/SILVA/SILVA16sDNA.taxa"
df_silvaSpecies = pd.read_csv(file, sep='\t', dtype=str, header=None)

def getSpecies(taxid):
    '''
    return species id with species id or subspecies id
    '''
    if taxid not in dc_rank:
        return None
    if dc_rank[taxid] == 'species':
        return taxid
    n = 0
    while n<5:
        n += 1
        taxid = dc_node[taxid]
        if dc_rank[taxid] == 'species':
            return taxid
    #print('something wrong')
    return None

df_silvaSpecies['taxid'] = df_silvaSpecies[1].apply(lambda x:getSpecies(x))
df_silvaSpecieskeep = df_silvaSpecies[df_silvaSpecies.taxid.notnull()]
df_silvaSpecieskeep['changed'] = df_silvaSpecieskeep.apply(lambda x:x[1]!=x['taxid'], axis=1)
df_silvaSpecieskeep.to_csv("/mnt/d/linux/W/NCBI/SILVA/SILVA16sDNA.taxaSubspeciesCorrected",columns=[0,'taxid'],header=None,index=None,sep='\t')

silva_taxids = set(df_silvaSpecies['taxid'])

print(len([e for e in taxids if e in silva_taxids]))
## read ezbio 16s
file = "/mnt/d/linux/W/NCBI/Ezbio/eztaxon_qiime_full.taxa"
df_ezbioSpecies = pd.read_csv(file, sep='\t', dtype=str, header=None)
ezbio_taxids = set(df_ezbioSpecies[1])
print(len([e for e in taxids if e in ezbio_taxids]))

print(len([e for e in taxids if e in ezbio_taxids or e in silva_taxids]))

## extract 16s sequences from ezbio and silva
taxid_with16s = [e for e in taxids if e in ezbio_taxids or e in silva_taxids]
taxid_no16s = [e for e in taxids if e not in taxid_with16s]
taxid_no16sInsilva = [e for e in taxids if e not in silva_taxids]

df_silvaSpecies['keep'] = df_silvaSpecies['taxid'].apply(lambda x:x in taxids_s)
df_silvaKeep = df_silvaSpecies[df_silvaSpecies['keep']]
silva_Keep = set(df_silvaKeep[0])

from collections import Counter
a = Counter(df_silvaKeep[1]).most_common()
file = '/mnt/d/linux/W/NCBI/SILVA/SILVA16sDNA.fasta'
fileout = '/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018SILVA.fasta'
fout = open(fileout,'w')
from Bio import SeqIO
for seq in SeqIO.parse(open(file),'fasta'):
    if seq.id.split('.')[0] in silva_Keep:
        fout.write('>'+seq.id+'\n'+str(seq.seq)+'\n')
fout.close()

## remove duplicate sequences in Human_gut_bacterial_species_2018SILVA.fasta
import sys
sys.path.append('/mnt/c/Users/ATPs/Documents/GitHub/XCProject/fasta/')
import largeFastaFile
fileout = '/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018SILVA.fasta'
ls_SILVA = largeFastaFile.open_fasta_to_list(fileout)
ls_SILVAuni = largeFastaFile.fasta_uni_keepone(ls_SILVA)
largeFastaFile.saveFastaListToFile(ls_SILVAuni,'/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018SILVAuni.fasta')
dc_id2tax = {e[1][0]:e[1]['taxid'] for e in df_silvaKeep.iterrows()}
dc_tax2seq = {}
file = '/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018SILVAuni.fasta'
for seq in SeqIO.parse(open(file),'fasta'):
    seq_acc = seq.id.split('.')[0]
    seq_tax = dc_id2tax[seq_acc]
    if seq_tax not in dc_tax2seq:
        dc_tax2seq[seq_tax] = []
    dc_tax2seq[seq_tax].append(seq)
dc_tax2sequni = {_k:largeFastaFile.fasta_within_seq_big_faster(_v) for _k, _v in dc_tax2seq.items()}

ls_SILVAnowithin = [j for i in dc_tax2sequni.values() for j in i]
largeFastaFile.saveFastaListToFile(ls_SILVAnowithin,'/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018SILVAbest.fasta')


df_ezbioKeep = df_ezbioSpecies[df_ezbioSpecies[1].apply(lambda x:x in taxids_s)]
b = Counter(df_ezbioKeep[1]).most_common()
ezbio_keep = set(df_ezbioKeep[0])
file = '/mnt/d/linux/W/NCBI/Ezbio/eztaxon_qiime_full.fasta'
fileout = '/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018Ezbio.fasta'
fout = open(fileout,'w')
from Bio import SeqIO
for seq in SeqIO.parse(open(file),'fasta'):
    if seq.id in ezbio_keep:
        fout.write('>'+seq.id+'\n'+str(seq.seq)+'\n')
fout.close()

## find ezbio seqid
txt = open('/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018SILVAbest.fasta').read()
ls_ez = largeFastaFile.open_fasta_to_list('/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018Ezbio.fasta')
ls_ezOnly = [e for e in ls_ez if str(e.seq) not in txt]
largeFastaFile.saveFastaListToFile(ls_ezOnly,'/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018EzbioOnly.fasta')

# blast ezOnly against SILVA, keep those with identity less than 99.7%. 39 seqs left
ls_ezFurtherCheck = '141783 141832 86987 97235 141769 80986 87793 141069 141068 142894 89066 88737 85921 141811 142911 143612 100887 88565 89143 141814 141125 141124 93691 85844 89067 141786 143344 89154 95759 80693 142574 142575 90380 140608 80815 114768 89954 91907 142555'.split()
dc_ez = largeFastaFile.open_fasta_to_dic('/mnt/d/linux/W/NCBI/Ezbio/eztaxon_qiime_full.fasta')
ls_ezFurtherCheckseq = [dc_ez[e] for e in ls_ezFurtherCheck]
largeFastaFile.saveFastaListToFile(ls_ezFurtherCheckseq,'/mnt/d/linux/W/NCBI/MoonRef/Human_gut_bacterial_species_2018EzbioFurtherCheck.fasta')

## searched NCBI refseq genomes, no identical sequences. Change those 39 sequences to NCBI seqs. 
## searched NCBI nr
## too complex. Use another way.

## check if all target species in nt_micro
file = "/mnt/d/linux/W/NCBI/nt_ref/nt_micro.taxa"
df_ntmicroSpecies = pd.read_csv(file, sep='\t', dtype=str, header=None)
ntmicro_taxids = set(df_ntmicroSpecies[1])
print(len([e for e in taxids if e in ntmicro_taxids]))

## check if all target species in speicesLineage
file = "/mnt/d/linux/W/NCBI/taxonamy/speicesLineage"
df_Lineage = pd.read_csv(file, sep='\t', dtype=str)
Linieage_taxids = set(df_Lineage["species"])
print(len([e for e in taxids if e in Linieage_taxids]))


