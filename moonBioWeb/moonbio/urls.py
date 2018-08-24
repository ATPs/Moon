"""moonbio URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from . import views


from django.conf.urls.static import static
from django.conf import settings

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.index, name='main page'),
    path('blast', views.blast, name='blast search'),
    path('align', views.align, name='align nucleotide or protein sequences. 比对两条（组）DNA或者蛋白质序列'),
    path('human_microbes', views.human_microbes, name='taxonamy for human_microbes'),
    path('env_microbes', views.env_microbes, name='taxonamy for env_microbes'),
    path('env_fungi', views.env_fungi, name='taxonamy for env_fungi'),
    path('this_month', views.this_month, name='result for this month'),
    path(r'\w*MOONINDEX', views.index, name='main page'),
    path('blastntaxa/', include('blastntaxa.urls')),
    path('helps/', include('helps.urls')),
]


urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)



