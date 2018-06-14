from django.urls import path

from . import views

urlpatterns = [
    path('', views.upload_file, name='index'),
    path('abi', views.abi_converter, name='abi'),
    path('today', views.today,name='today'),
    path('contact/',views.contact, name = 'contact'), 
    path('index.html',views.model_form_upload, name='upload')
]

