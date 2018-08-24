from django.shortcuts import render
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from .forms import UploadFileForm, DocumentForm
from . import blastntaxa
import datetime
import os

handle_uploaded_file = blastntaxa.handle_uploaded_file

# Create your views here.
def index(request):
    return render(request, 'blastntaxa/home.html')

def index2(request):
    return render(request, 'index.html')

def blastn(request):
    if request.method == 'POST' and request.FILES['myfile']:
        myfile = request.FILES['myfile']
        fs = FileSystemStorage()
        filename = fs.save(myfile.name, myfile)
        uploaded_file_url = fs.url(filename)
        return render(request, 'static/blastn.html', {
            'uploaded_file_url': uploaded_file_url
        })
    return render(request, 'static/blastn.html')

def contact(request):
    return render(request, 'blastntaxa/basic.html', {'content':["If you would like to contact me, please email me","caoxl@moonbio.com"]})
    

def upload_file(request):
    if request.method == 'POST':
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            fs = FileSystemStorage()
            filename = fs.save(myfile.name, myfile)
            filename = myfile.name
        a = datetime.datetime.now()
        query_pre = "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
        query_name = '/mnt/d/linux/M/www/Django/moonbio/media/' + query_pre
        sequences = request.POST.get('sequences')
        #print(sequences)
        if len(sequences)>2:
            open(query_name,'w').write(sequences)
            filename = query_pre
        databases = request.POST.getlist('database')
        taxa_order = request.POST.get('taxonomy')
        #print(databases, myfile,taxa_order)
        #return render(request, 'blastntaxa/home.html')
        if len(databases) == 0:
            return render(request, 'blastntaxa/home.html', {'result':'failed! No database selected!'})
        foldername, best_table = handle_uploaded_file(filename, databases = databases,taxa_order=taxa_order)
        #print(best_table)
        return render(request, 'blastntaxa/home.html', {
            'uploaded_file_url': '/media/'+foldername+'.zip','best_table':best_table
        })
    return render(request, 'blastntaxa/home.html')

def today(request):
    
    a = datetime.datetime.now()
    foldernameYearMonth = "%0d%02d"%(a.year,a.month)
    path = os.path.join('/mnt/d/linux/M/www/Django/moonbio/media/','taxonomy',foldernameYearMonth)
    files = os.listdir(path)
    a = datetime.datetime.now()
    this_month = "%0d%02d"%(a.year,a.month)
    files = [os.path.join('taxonomy',foldernameYearMonth,e) for e in files if this_month in e]
    #print(files)
    return render(request, 'blastntaxa/today.html', {'content':files})

def model_form_upload(request):
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            return redirect('home')
    else:
        form = DocumentForm()
    return render(request, 'blastntaxa/model_form_upload.html', {
        'form': form
    })

def abi_converter(request):
    if request.method == 'POST' and request.FILES['myfile']:
        myfile = request.FILES['myfile']
        filename = myfile.name
        fs = FileSystemStorage()
        fs.save(myfile.name, myfile)
        merge = request.POST.get('merge')
        if merge == '0':
            merge = False
        else:
            merge = True
        if len(filename) < 5:
            return render(request, 'blastntaxa/abi_converter.html', {'result':'failed! file name is wrong!'})
        if filename[-4:] not in {'.ab1','.AB1','.zip','.ZIP'}:
            return render(request, 'blastntaxa/abi_converter.html', {'result':'failed! file extension is not .zip or .ab1'})
        if filename[-4:] == '.ab1' or filename[-4:] == '.AB1':
            result, seq = blastntaxa.abi_converter_single(filename)
            return render(request, 'blastntaxa/abi_converter.html', {'result':result,'fasta':seq})
        if filename[-4:] == '.zip' or filename[-4:] == '.ZIP':
            if not merge:
                result, seq, finalpath = blastntaxa.abi_converter(filename)
                result += '\nResults:'
                return render(request, 'blastntaxa/abi_converter.html', {'result':result,'fasta':seq,'uploaded_file_url': finalpath})
            else:
                result, seq, finalpath = blastntaxa.abi_converter(filename)
                result2, seq = blastntaxa.seqmerger(finalpath)
                result = result + '\nAfter merge, some low quality ones:\n'+result2
                result += '\nResults:'
                return render(request, 'blastntaxa/abi_converter.html', {'result':result,'fasta':seq,'uploaded_file_url': finalpath})
    return render(request, 'blastntaxa/abi_converter.html')