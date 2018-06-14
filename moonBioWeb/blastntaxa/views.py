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
    if request.method == 'POST' and request.FILES['myfile']:
        myfile = request.FILES['myfile']
        databases = request.POST.getlist('database')
        taxa_order = request.POST.get('taxonomy')
        #print(databases, myfile,taxa_order)
        #return render(request, 'blastntaxa/home.html')
        if len(databases) == 0:
            return render(request, 'blastntaxa/home.html', {'result':'failed! No database selected!'})
        fs = FileSystemStorage()
        filename = fs.save(myfile.name, myfile)
        foldername = handle_uploaded_file(myfile.name, databases = databases,taxa_order=taxa_order)
        return render(request, 'blastntaxa/home.html', {
            'uploaded_file_url': '/media/'+foldername+'.zip'
        })
    return render(request, 'blastntaxa/home.html')

def today(request):
    path = '/mnt/d/linux/M/www/Django/moonbio/media/'
    files = os.listdir(path)
    a = datetime.datetime.now()
    this_month = "%0d%02d"%(a.year,a.month)
    files = [e for e in files if this_month in e]
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
        if len(filename) < 5:
            return render(request, 'blastntaxa/abi_converter.html', {'result':'failed! file name is wrong!'})
        if filename[-4:] not in {'.ab1','.AB1','.zip','.ZIP'}:
            return render(request, 'blastntaxa/abi_converter.html', {'result':'failed! file extension is not .zip or .ab1'})
        if filename[-4:] == '.ab1' or filename[-4:] == '.AB1':
            result, seq = blastntaxa.abi_converter_single(filename)
            return render(request, 'blastntaxa/abi_converter.html', {'result':result,'fasta':seq})
        if filename[-4:] == '.zip' or filename[-4:] == '.ZIP':
            result, seq = blastntaxa.abi_converter(filename)
            result += '\nResults:'
            return render(request, 'blastntaxa/abi_converter.html', {'result':result,'fasta':seq,'uploaded_file_url': '/media/'+filename[:-4]+'_fas'+'.zip'})
    return render(request, 'blastntaxa/abi_converter.html')