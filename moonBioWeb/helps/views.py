from django.shortcuts import render

# Create your views here.
def index(request):
    return render(request, 'helps/helps.html')

def viewTest(request):
    return render(request, 'helps/test.html')