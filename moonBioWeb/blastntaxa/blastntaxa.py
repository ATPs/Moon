import datetime
import os
from django.core.files.storage import FileSystemStorage
import zipfile
import shutil
import subprocess
import sys
import pandas as pd
sys.path.append('/mnt/d/linux/M/www/Django/TaxaFinder/')

def handle_uploaded_file(filename, databases, taxa_order = '0'):
    a = datetime.datetime.now()
    foldernameYearMonth = "%0d%02d"%(a.year,a.month)
    foldername = "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
    
    folder = os.path.join('/mnt/d/linux/W/MoonNt/', foldernameYearMonth, foldername)
    if not os.path.exists(folder):
        os.makedirs(folder)
    #os.chdir(folder)
    os.rename(os.path.join('/mnt/d/linux/M/www/Django/moonbio/media',filename), os.path.join(folder, filename))
    if '.zip' in filename or '.ZIP' in filename:
        filezip = zipfile.ZipFile(os.path.join(folder, filename))
        filezip.extractall(os.path.join(folder, 'seqs'))
    else:
        if not os.path.exists(os.path.join(folder, 'seqs')):
            os.mkdir(os.path.join(folder, 'seqs'))
        shutil.move(os.path.join(folder, filename), os.path.join(folder, 'seqs',filename+'.seq'))
    #os.remove(filename)
    best_table = '<h5>我只想看看最佳结果：</h5>'
    for db in databases:
        commandline = '/home/p/anaconda/anaconda3_5.2.0/bin/python3.6 /mnt/d/linux/M/www/Django/moonbio/blastntaxa/TaxaFinder_Moon20180814.py {folder} --db {db} --taxa_order {taxa_order}'.format(folder=folder, db=db, taxa_order=taxa_order)
        subprocess.run(commandline,shell=True)
        file_best = os.path.join(folder, foldername+'_Best'+db+'Result.xlsx')
        _df = pd.read_excel(file_best)
        best_table += '<p>{db}</p>'.format(db=db)
        best_table += _df.to_html(columns=['query', 'subject', 'length', 'matchLength', 'identical', 'qcover', 'identity', 'taxID', 'species', 'genus', 'family', 'order', 'class', 'phylum'],index=False,justify='center')
    #shutil.rmtree(folder+'/seqs/')
    html_folder = os.path.join('/mnt/d/linux/M/www/Django/moonbio/media','taxonomy',foldernameYearMonth)
    if not os.path.exists(html_folder):
        os.makedirs(html_folder)
    if filename.endswith('.zip'):
        filename = filename[:-4]
    shutil.make_archive(os.path.join(html_folder,foldername+ '_'.join(databases)+'_'+filename), 'zip', folder)
    return os.path.join('taxonomy',foldernameYearMonth,foldername+ '_'.join(databases)+'_'+filename), best_table

def abi_converter_single(filename):
    import abi
    filename = os.path.join('/mnt/d/linux/M/www/Django/moonbio/media',filename)
    seq =abi.cleanABI(filename, qmin = 15, window = 10, minlen = 200, outputfile = False)
    os.remove(filename)
    if seq == '':
        return "quality low. No output", ''
    return 'Below is the suquence:', seq

def abi_converter(filename):
    import abi
    import glob
    filepath = os.path.join('/mnt/d/linux/M/www/Django/moonbio/media',filename)
    a = datetime.datetime.now()
    foldername = "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
    foldernameYearMonth = "%0d%02d"%(a.year,a.month)
    folder = os.path.join('/mnt/d/linux/W/MoonNt/AB1/', foldername)
    if not os.path.exists(folder):
        os.makedirs(folder)
    #os.chdir(folder)
    os.rename(filepath, os.path.join(folder, filename))
    filezip = zipfile.ZipFile(os.path.join(folder, filename))
    filezip.extractall(folder)
    files = glob.glob(folder+"/**/*.ab1",  recursive=True)
    files += glob.glob(folder+"/**/*.AB1",  recursive=True)
    seqs = []
    nogood = ['low quality files:']
    for f in files:
        seq =abi.cleanABI(f, qmin = 15, window = 10, minlen = 200, outputfile = False)
        if seq == '':
            nogood.append(os.path.basename(f))
        else:
            seqs.append(seq)
        os.remove(f)
        #print(seq)
    
    htmlpath = os.path.join('/mnt/d/linux/M/www/Django/moonbio/media', 'AB1',foldernameYearMonth)
    os.rename(os.path.join(folder, filename),os.path.join(htmlpath, foldername + filename))
    if not os.path.exists(htmlpath):
        os.makedirs(htmlpath)
    fout = open(os.path.join(htmlpath, foldername + filename[:-4]+'.seq'),'w')
    [fout.write(e) for e in seqs]
    fout.close()
    shutil.rmtree(folder)
    finalpath = os.path.join('/media/AB1',foldernameYearMonth,foldername + filename[:-4]+'.seq')
    return ' '.join(nogood), ''.join(seqs), finalpath



def seqmerger(filename):
    import abi
    filename1 = '/mnt/d/linux/M/www/Django/moonbio'+filename
    print(filename1)
    return abi.mergerFastaSeqInFile(filename1)