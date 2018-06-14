import datetime
import os
from django.core.files.storage import FileSystemStorage
import zipfile
import shutil
import subprocess
import sys
sys.path.append('/mnt/c/Users/ATPs/Documents/GitHub/Moon/TaxaFinder/')

def handle_uploaded_file(filename, databases, taxa_order = '0'):
    a = datetime.datetime.now()
    foldername = "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
    folder = '/mnt/d/linux/W/MoonNt/django/' + foldername
    if not os.path.exists(folder):
        os.makedirs(folder)
    os.chdir(folder)
    os.rename(os.path.join('/mnt/d/linux/M/www/Django/moonbio/media',filename), os.path.join(folder, filename))
    if '.zip' in filename or '.ZIP' in filename:
        filezip = zipfile.ZipFile(filename)
        filezip.extractall('./seqs/')
    else:
        if not os.path.exists('seqs'):
            os.mkdir('seqs')
        shutil.copy(filename, './seqs/'+filename+'.seq')
    for db in databases:
        commandline = '/mnt/d/linux/P/anaconda3/bin/python3.6 /mnt/c/Users/ATPs/Documents/GitHub/Moon/TaxaFinder_Moon20180330 {folder} --db {db} --taxa_order {taxa_order}'.format(folder=folder, db=db, taxa_order=taxa_order)
        subprocess.run(commandline,shell=True)
    shutil.make_archive(os.path.join('/mnt/d/linux/M/www/Django/moonbio/media',foldername+ '_'.join(databases)), 'zip', folder)
    return foldername+ '_'.join(databases)

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
    folder = '/mnt/d/linux/W/MoonNt/django/' + foldername
    if not os.path.exists(folder):
        os.makedirs(folder)
    os.chdir(folder)
    os.rename(filepath, os.path.join(folder, filename))
    filezip = zipfile.ZipFile(filename)
    filezip.extractall('./')
    files = glob.glob(folder+"/**/*.ab1",  recursive=True)
    files += glob.glob(folder+"/**/*.AB1",  recursive=True)
    seqs = []
    nogood = ['low quality files:']
    for f in files:
        seq =abi.cleanABI(f, qmin = 15, window = 10, minlen = 200, outputfile = True)
        if seq == '':
            nogood.append(os.path.basename(f))
        else:
            seqs.append(seq)
        os.remove(f)
        #print(seq)
    os.remove(filename)
    shutil.make_archive(filepath[:-4]+'_fas', 'zip', folder)
    shutil.rmtree(folder)
    return ' '.join(nogood), ''.join(seqs)
