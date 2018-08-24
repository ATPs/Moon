from django.shortcuts import render
from django.shortcuts import render_to_response
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponse
from django.conf import settings
import os
import datetime
import sys

# Create your views here.



def index(request):
    #return render_to_response('static/index.html')
    #return render(request, 'static/index.html')
    #html = open(os.path.join(settings.BASE_DIR,'static/index.html')).read()
    return render(request, 'index.html')



def blast(request):
    if request.method == 'POST':
        a = datetime.datetime.now()
        query_pre = "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
        query_name = '/mnt/d/linux/M/www/Django/moonbio/media/blast/' + query_pre
        #print('request.FILES',request.FILES)
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            if myfile:
                #print('myfile is',myfile)
                fs = FileSystemStorage()
                fs.save(query_name, myfile)
        sequences = request.POST.get('sequences')
        #print(sequences)
        if len(sequences)>2:
            open(query_name,'w').write(sequences)
        
        blast_method = request.POST.get('blast_method')
        blast_db = request.POST.get('blast_db')
        outfmt = request.POST.get('blast_out')
        dcDB = {"fungi":"/mnt/d/linux/W/NCBI/fungi/ITS7.2",
        "SILVA16s":"/mnt/d/linux/W/NCBI/SILVA/SILVA16sDNA",
        "SILVA18s":"/mnt/d/linux/W/NCBI/SILVA/SILVA18sDNA",
        "Ezbio":"/mnt/d/linux/W/NCBI/Ezbio20180810/Ezbio20180810",
        #"nt_micro":"/mnt/d/linux/W/NCBI/nt_ref/nt_micro",
        "nt_micro":"/mnt/d/linux/W/NCBI/nt20180810/nt_micro",
        "nt":"/mnt/d/linux/W/NCBI/nt201803/nt",
        "nr":"/mnt/d/linux/W/NCBI/nr201804/nr"
        }
        blast_db = dcDB[blast_db]
        outname = query_name + '.html'
        #print('blast_method', blast_method,'blast_db', blast_db)
        commandline = '/home/p/blast/bin/{blast_method} -query {query_name} -db {blast_db} -num_threads 24  -outfmt {outfmt} -out {outname} -html -num_alignments 100  '.format(blast_method=blast_method, blast_db=blast_db, query_name=query_name, outname=outname, outfmt=outfmt)
        import subprocess
        subprocess.run(commandline,shell=True)
        blast_result = open(outname).read()
        os.remove(outname)
        os.remove(query_name)
        if outfmt == '0':
            return HttpResponse(blast_result)
        else:
            return HttpResponse(blast_result,content_type='text/plain; charset=utf8')
    return render(request, 'blast.html')

def align(request):
    if request.method == 'POST':
        folder = '/mnt/d/linux/M/www/Django/moonbio/media/blast/'
        a = datetime.datetime.now()
        query_name = folder + "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)+"query"
        subject_name = folder + "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)+"subject"
        outname = folder + "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)+"out"
        
        #print('request.FILES',request.FILES)
        if 'myfile1' in request.FILES:
            myfile1 = request.FILES['myfile1']
            if myfile1:
#                print('myfile1 is',myfile1)
                fs = FileSystemStorage()
                fs.save(query_name, myfile1)
        sequences1 = request.POST.get('sequences1')
#        print(sequences1)
        if len(sequences1)>2:
            open(query_name,'w').write(sequences1)
        
        if 'myfile2' in request.FILES:
            myfile2 = request.FILES['myfile2']
            if myfile2:
#                print('myfile2 is',myfile2)
                fs = FileSystemStorage()
                fs.save(subject_name, myfile2)
        sequences2 = request.POST.get('sequences2')
#        print(sequences2)
        if len(sequences2)>2:
            open(subject_name,'w').write(sequences2)
        
        blast_method = request.POST.get('blast_method')
        
        commandline = '/home/p/blast/bin/{blast_method} -query {query_name} -subject {subject_name}  -out {outname} -html'.format(blast_method=blast_method, subject_name=subject_name, query_name=query_name, outname=outname)
        import subprocess
        subprocess.run(commandline,shell=True)
        blast_result = open(outname).read()
        os.remove(outname)
        os.remove(query_name)
        os.remove(subject_name)
        return HttpResponse(blast_result)
    return render(request, 'align.html')

def human_microbes(request):
    if request.method == 'POST':
        folder = '/mnt/d/linux/M/www/Django/moonbio/media/HumanMicrobes'
        outexcelname = 'MNH_human_micro_all'
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            fs = FileSystemStorage()
            filenameStore = os.path.join(folder,'input/',myfile.name)
            filename = fs.save(filenameStore, myfile)
            #filename = myfile.name
        
        print(filename, filenameStore)
        os.rename(filename,filenameStore)
        
        sys.path.append('/mnt/d/linux/M/www/Django/TaxaFinder/')
        import HumanMicrobesTaxaFinder_Excel_Moon20180814
        
        foldernameYearMonthDay = HumanMicrobesTaxaFinder_Excel_Moon20180814.taxaFinder(filenameStore,folder = folder,outexcelname = outexcelname)
        return render(request, 'blastntaxa/bacteria_taxa.html', {
            'uploaded_file_url': '/media/HumanMicrobes/'+foldernameYearMonthDay+outexcelname+'.xlsx',
        })
    return render(request, 'blastntaxa/bacteria_taxa.html')

def env_microbes(request):
    if request.method == 'POST':
        folder = '/mnt/d/linux/M/www/Django/moonbio/media/EnvMicrobes'
        outexcelname = 'MN_env_micro_all'
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            fs = FileSystemStorage()
            filenameStore = os.path.join(folder,'input/',myfile.name)
            filename = fs.save(filenameStore, myfile)
            #filename = myfile.name
        
        print(filename, filenameStore)
        os.rename(filename,filenameStore)
        
        sys.path.append('/mnt/d/linux/M/www/Django/TaxaFinder/')
        import HumanMicrobesTaxaFinder_Excel_Moon20180814
        foldernameYearMonthDay = HumanMicrobesTaxaFinder_Excel_Moon20180814.taxaFinder(filenameStore,folder = folder,outexcelname = outexcelname)
        return render(request, 'blastntaxa/env_bacteria_taxa.html', {
            'uploaded_file_url': '/media/EnvMicrobes/'+foldernameYearMonthDay+outexcelname+'.xlsx',
        })
    return render(request, 'blastntaxa/env_bacteria_taxa.html')

def env_fungi(request):
    if request.method == 'POST':
        folder = '/mnt/d/linux/M/www/Django/moonbio/media/EnvFungi'
        outexcelname = 'MN_env_fungi_all'
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            fs = FileSystemStorage()
            filenameStore = os.path.join(folder,'input/',myfile.name)
            filename = fs.save(filenameStore, myfile)
            #filename = myfile.name
        
        print(filename, filenameStore)
        os.rename(filename,filenameStore)
        
        sys.path.append('/mnt/d/linux/M/www/Django/TaxaFinder/')
        import HumanMicrobesTaxaFinder_Excel_Moon20180814
        foldernameYearMonthDay = HumanMicrobesTaxaFinder_Excel_Moon20180814.taxaFinder(filenameStore,folder = folder,outexcelname = 'MN_env_fungi_all', dbs=['fungi','nt_micro'])
        return render(request, 'blastntaxa/env_fungi_taxa.html', {
            'uploaded_file_url': '/media/EnvFungi/'+foldernameYearMonthDay+outexcelname+'.xlsx',
        })
    return render(request, 'blastntaxa/env_fungi_taxa.html')
 
def this_month(request):
    
    a = datetime.datetime.now()
    foldernameYearMonth = "%0d%02d"%(a.year,a.month)
    path = os.path.join('/mnt/d/linux/M/www/Django/moonbio/media/')
    import glob
    files = glob.glob(path+'*/*.xlsx')
    files = [e.split('www/Django/moonbio/media/')[1] for e in files if foldernameYearMonth in e]
    #print(files)
    return render(request, 'blastntaxa/today.html', {'content':files})