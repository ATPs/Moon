
# coding: utf-8

# # Pipeline to Assemble the genome of bacteria
# 
# Combination of two bacteria genome assembling pipelines, Tychus and Mypro.  
# Realization of Tychus with Python scripts. Add SOAPdenovo, improve the genus database for prokka with data from Mypro https://sourceforge.net/projects/sb2nhri/files/MyPro/ 

# ## setup the files and parameters you are going to work with.

# ### two reads files

# In[1]:

import sys
read1 = sys.argv[1]
read2 = sys.argv[2]


# ### parameters

# In[54]:


thread = 24
outfolder = sys.argv[3]
bacteriaStrain = sys.argv[4] # change as needed
cleanOrNot = True # whether or not to remove intermediate files


# ### usually fixed parameters

# In[3]:


adapters = '/home/p/Trimmomatic-0.38/adapters/TruSeq3-PE.fa'


# ## import required packages

# In[53]:


import subprocess
import os
import shutil
from Bio import SeqIO
import time


# ## set up some constants
# Below are some programs that will be used

# In[5]:


## program fastqc
FASTQC = '/home/p/FastQC/fastqc'
## program, Trim reads
TRIMMOMATIC = '/home/p/Trimmomatic-0.38/trimmomatic-0.38.jar'
## program, reads quality report and Trim
FASTP = '/home/p/fastp/fastp'
## program, Identify the best kmer for genome assembling
KMERGENIE = '/home/p/KmerGenie/kmergenie-1.7048/kmergenie'
## program, assemble the genome
ABYSS = '/home/p/anaconda/anaconda3_5.2.0/bin/abyss-pe'
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
CISAPATH = '/home/p/CISA/CISA1.3/'
CISA = '/home/p/CISA/CISA1.3/CISA.py'
## program, merge assembled results
NOVO_STITCH = '/home/p/Novo_Stitch/scripts/main.py'
REFALIGNER = '/home/p/RefalignerAssembler/latest/RefAligner'
GLPK = '/usr/local/bin/glpsol'
## program, annotate merged genome
PROKKA = '/home/p/prokka/prokka-master/bin/prokka'
## program, evaluate result
QUAST = '/home/p/quast/quast-4.6.3/quast.py'

## used programs
NCUMER = '/home/p/quast/quast-4.6.3/quast_libs/MUMmer/nucmer'
MAKEBLASTDB = '/home/p/blast/bin/makeblastdb'
BLASTN = '/home/p/blast/bin/blastn'


# ## check quality of input reads with fastqc

# creat the folder to store fastqc result

# In[6]:


if not os.path.exists(outfolder):
    os.mkdir(outfolder)
os.chdir(outfolder)
if not os.path.exists('fastqc'):
    os.mkdir('fastqc')


# run the command

# In[7]:


command_line = '{FASTQC} {read1} {read2} -t {thread} -o {outfolder}/fastqc'.format(                FASTQC = FASTQC, read1 = read1, read2 = read2, thread = thread, outfolder = outfolder)
print(command_line)
subprocess.run(command_line,shell=True)


# ## trim reads with TRIMMOMATIC

# creat the folder to store the result

# In[8]:


if not os.path.exists('Trimmomatic'):
    os.mkdir('Trimmomatic')


# In[9]:


command_line = '''java -jar {TRIMMOMATIC} PE -threads {thread} -phred33 {read1} {read2}                 ./Trimmomatic/read1.pair.fq ./Trimmomatic/read1.single.fq ./Trimmomatic/read2.pair.fq ./Trimmomatic/read2.single.fq                 ILLUMINACLIP:{adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 '''.format(                TRIMMOMATIC = TRIMMOMATIC, read1 = read1, read2 = read2, thread = thread, adapters = adapters)
print(command_line)
result = subprocess.run(command_line,shell=True, check=True)
print(result)


# ## Identify best kmer with kmergenie

# In[10]:


if not os.path.exists('kmergenie'):
    os.mkdir('kmergenie')
fout = open('readlist.txt','w')
fout.write(os.path.join(outfolder,'Trimmomatic','read1.pair.fq') + '\n')
fout.write(os.path.join(outfolder,'Trimmomatic','read2.pair.fq') + '\n')
fout.close()


# In[11]:


command_line = '''{KMERGENIE} {outfolder}/readlist.txt -t {thread} -o {outfolder}/kmergenie/kmergenie.histogram '''.format(                KMERGENIE = KMERGENIE, outfolder = outfolder, thread = thread, adapters = adapters)
print(command_line)
result = subprocess.check_output(command_line, shell=True)


# In[12]:


result = result.decode('utf-8')
print(result)


# In[13]:


result.split('\n')[-5]
bestkmer = result.split('\n')[-5].split()[-1]
print('best kmer is ', bestkmer)


# ## run abyss

# In[14]:


if not os.path.exists('abyss'):
    os.mkdir('abyss')
os.chdir('abyss')


# In[15]:


command_line = '''{ABYSS} k={bestkmer} np={thread}  in="{read1} {read2}" name=abyss'''.format(                ABYSS = ABYSS, bestkmer=bestkmer, outfolder = outfolder, thread = thread,                read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq'))
print(command_line)
result = subprocess.check_output(command_line, shell=True)


# In[16]:


result = result.decode('utf-8')
print(result)
os.chdir('..')


# ## run velvet

# In[17]:


if not os.path.exists('velvet'):
    os.mkdir('velvet')
os.chdir('velvet')


# In[18]:


os.getcwd()


# In[19]:


command_line = '''{VELVETH} . {bestkmer} -separate -fastq -shortPaired {read1} {read2} '''.format(                VELVETH = VELVETH, bestkmer=bestkmer, outfolder = outfolder, thread = thread,                read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq'))
print(command_line)
result = subprocess.check_output(command_line, shell=True)


# In[20]:


result = result.decode('utf-8')
print(result)


# In[21]:


command_line = '''{VELVETG} .  -exp_cov auto -cov_cutoff auto '''.format( VELVETG = VELVETG)
print(command_line)
result = subprocess.check_output(command_line, shell=True)


# In[22]:


result = result.decode('utf-8')
print(result)
os.chdir('..')


# ## run SPAdes

# In[23]:


if not os.path.exists('spades'):
    os.mkdir('spades')
os.chdir('spades')


# In[24]:


os.getcwd()


# In[25]:


command_line = '''{SPADES} --pe1-1 {read1} --pe1-2  {read2} -t {thread} -o .  '''.format(                SPADES = SPADES, bestkmer=bestkmer, outfolder = outfolder, thread = thread,                read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq'))
print(command_line)
result = subprocess.check_output(command_line, shell=True)


# In[26]:


result = result.decode('utf-8')
print(result)
os.chdir('..')


# ## run IDBA-UD

# In[27]:


if not os.path.exists('IDBA-UD'):
    os.mkdir('IDBA-UD')
os.chdir('IDBA-UD')


# In[28]:


command_line = '''{FQ2FA} --merge --filter {read1} {read2} idba-paired-contigs.fa '''.format(                FQ2FA = FQ2FA, bestkmer=bestkmer, outfolder = outfolder, thread = thread,                read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq'))
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
print(result)


# In[29]:


command_line = '''{IDBA_UD} -r idba-paired-contigs.fa --num_threads {thread} -o . '''.format(                IDBA_UD = IDBA_UD, bestkmer=bestkmer, outfolder = outfolder, thread = thread,                read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq'))
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
print(result)


# In[30]:


os.chdir('..')


# ## run soapdenovo

# In[31]:


if not os.path.exists('soapdenovo'):
    os.mkdir('soapdenovo')
os.chdir('soapdenovo')


# In[32]:


fw=open('myconfig','w')
fw.write('max_rd_len=150\n[LIB]\navg_ins=350\nreverse_seq=0\nasm_flags=3\nq1={read1}\nq2={read2}\n'.format(        read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq')))
fw.close()


# In[33]:


command_line = '''{SOAPDENOVO} all -s myconfig -K {bestkmer} -R -o soapdenovo -d 10 -p 24  1>ass.log 2>ass.err '''.format(                SOAPDENOVO = SOAPDENOVO, bestkmer=bestkmer, outfolder = outfolder, thread = thread,                read1 = os.path.join(outfolder,'Trimmomatic','read1.pair.fq'), read2 = os.path.join(outfolder,'Trimmomatic','read2.pair.fq'))
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
print(result)


# In[34]:


os.chdir('..')


# ## run CISA

# In[35]:


if not os.path.exists('CISA'):
    os.mkdir('CISA')
os.chdir('CISA')


# In[36]:


shutil.copy('../abyss/abyss-contigs.fa','abyss.fa')
shutil.copy('../velvet/contigs.fa','velvet.fa')
shutil.copy('../spades/contigs.fasta','spades.fa')
shutil.copy('../IDBA-UD/contig.fa','IDBA-UD.fa')
shutil.copy('../soapdenovo/soapdenovo.contig','soapdenovo.fa')


# In[37]:


ctgs = ['soapdenovo.fa', 'velvet.fa', 'abyss.fa', 'spades.fa', 'IDBA-UD.fa']


# In[38]:


fw=open('ToMerge.config','w')
fw.write('count=5\n')
for i in ctgs:
    fw.write('data={fastafile},title={fastaname}\n'.format(fastafile=i,fastaname=i.split('.')[0]))
fw.write('min_length=400\nMaster_file=Merged.fa\nGap=1\n')
fw.close()


# In[39]:


command_line = '''{MERGE} ToMerge.config '''.format(MERGE = MERGE)
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
print(result)


# In[40]:


lis = open('Merge_info').readlines()
genomelens = [int(e.split()[-1]) for e in lis if 'Whole Genome' in e]
genomelen = max(genomelens)
print('max genome length is', genomelen)


# In[41]:


fw=open('CISA.config','w')
fw.write('''
genome={genomelen}
infile=Merged.fa
outfile=Final.fa
nucmer={NCUMER}
R2_Gap=0.95
CISA={CISAPATH}
makeblastdb={MAKEBLASTDB}
blastn={BLASTN}
'''.format(genomelen=genomelen,NCUMER=NCUMER,CISAPATH=CISAPATH,MAKEBLASTDB=MAKEBLASTDB,BLASTN=BLASTN))
fw.close()


# In[42]:


command_line = '''python {CISA} CISA.config '''.format(CISA = CISA)
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
#print(result)


# In[43]:


os.chdir('..')


# ## run prokka

# In[44]:


os.getcwd()


# In[45]:


if not os.path.exists('prokka'):
    os.mkdir('prokka')
os.chdir('prokka')


# In[46]:


ls_contigs = list(SeqIO.parse(open('../CISA/Final.fa'), 'fasta'))
outfile = '{bacteriaStrain}.genome.fa'.format(bacteriaStrain=bacteriaStrain)
fout = open(outfile,'w')
for _n, _e in enumerate(ls_contigs):
    fout.write('>'+bacteriaStrain+'.g'+str(_n+1)+'\n'+str(_e.seq)+'\n')
fout.close()


# In[47]:


command_line = '''{PROKKA} --kingdom Bacteria --cpus {thread}  {bacteriaStrain}.genome.fa --outdir .                 --prefix {bacteriaStrain} --force --metagenome --locustag {bacteriaStrain}'''.format(                PROKKA = PROKKA, thread=thread, bacteriaStrain=bacteriaStrain)
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
print(result)


# ## run quast

# In[48]:


os.chdir('..')


# In[50]:


os.getcwd()


# In[51]:


command_line = '''python {QUAST} -t {thread} -o quast                 ./prokka/{bacteriaStrain}.genome.fa ./CISA/abyss.fa.p.fa ./CISA/velvet.fa.p.fa ./CISA/spades.fa.p.fa                 ./CISA/IDBA-UD.fa.p.fa ./CISA/soapdenovo.fa.p.fa'''.format(                QUAST = QUAST, thread=thread, bacteriaStrain=bacteriaStrain)
print(command_line)
result = subprocess.check_output(command_line, shell=True)
result = result.decode('utf-8')
print(result)


# ## clean files

# In[56]:


if cleanOrNot:
    os.system('rm -rf {outfolder}/abyss'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/IDBA-UD'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/kmergenie'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/soapdenovo'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/spades'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/Trimmomatic'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/velvet'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/CISA'.format(outfolder=outfolder))
    os.system('rm -rf {outfolder}/readlist.txt'.format(outfolder=outfolder))
    

