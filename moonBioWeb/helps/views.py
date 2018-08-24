from django.shortcuts import render
from django.shortcuts import render_to_response
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponse
from django.conf import settings
import os
import datetime

# Create your views here.
def index(request):
    return render(request, 'helps/helps.html')

def viewTest(request):
    if request.method == 'POST':
        a = datetime.datetime.now()
        query_pre = "%0d%02d%02d%02d%02d%02d"%(a.year,a.month,a.day,a.hour,a.minute,a.second)
        
        #print('request.FILES',request.FILES)
        if 'myfile' in request.FILES:
            myfile = request.FILES['myfile']
            if myfile:
                print('myfile is',myfile)
                query_name = '/mnt/d/linux/M/www/Django/moonbio/media/blast/' + query_pre + myfile.name
                fs = FileSystemStorage()
                fs.save(query_name, myfile)
                query_name = '/mnt/d/linux/M/www/Django/moonbio/media/blast/' + query_pre + myfile.name
        sequences = request.POST.get('sequences')
        print(sequences)
        if len(sequences)>2:
            query_name = './media/blast/' + query_pre
            open(query_name,'w').write(sequences)
    return render(request, 'helps/test.html')