from django.shortcuts import render
from django.shortcuts import render_to_response
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponse
from django.conf import settings
import os
import datetime

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
        
        #print('request.FILES',request.FILES)
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            if myfile:
                #print('myfile is',myfile)
                query_name = './blast/' + query_pre + myfile.name
                fs = FileSystemStorage()
                fs.save(query_name, myfile)
                query_name = './media/blast/' + query_pre + myfile.name
        sequences = request.POST.get('sequences')
        #print(sequences)
        if len(sequences)>2:
            query_name = './media/blast/' + query_pre
            open(query_name,'w').write(sequences)
        
        blast_method = request.POST.get('blast_method')
        blast_db = request.POST.get('blast_db')
        dcDB = {"fungi":"/mnt/d/linux/W/NCBI/fungi/ITS7.2",
        "SILVA16s":"/mnt/d/linux/W/NCBI/SILVA/SILVA16sDNA",
        "SILVA18s":"/mnt/d/linux/W/NCBI/SILVA/SILVA18sDNA",
        "Ezbio":"/mnt/d/linux/W/NCBI/Ezbio/eztaxon_qiime_full",
        #"nt_micro":"/mnt/d/linux/W/NCBI/nt_ref/nt_micro",
        "nt_micro":"/mnt/c/P/DB/nt_micro/nt_micro",
        "nt":"/mnt/d/linux/W/NCBI/nt201803/nt",
        "nr":"/mnt/d/linux/W/NCBI/nr201804/nr"
        }
        blast_db = dcDB[blast_db]
        outname = query_name + '.html'
        #print('blast_method', blast_method,'blast_db', blast_db)
        commandline = '/mnt/d/linux/P/blast/bin/{blast_method} -query {query_name} -db {blast_db} -num_threads 8 -max_target_seqs 100 -out {outname} -html -num_alignments 100 '.format(blast_method=blast_method, blast_db=blast_db, query_name=query_name, outname=outname)
        import subprocess
        subprocess.run(commandline,shell=True)
        blast_result = open(outname).read()
        os.remove(outname)
        os.remove(query_name)
        return HttpResponse(blast_result)

    
    
    return render(request, 'blast.html')