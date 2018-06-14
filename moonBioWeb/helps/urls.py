from django.conf.urls import url, include
from . import views
from django.views.generic import ListView, DetailView
from helps.models import Post



urlpatterns = [
    #url(r'^$', views.index, name='index'),
    #url(r'$', ListView.as_view(queryset=Post.objects.all().order_by("-date"), template_name='helps/list.html')),
    url(r'^post$', ListView.as_view(queryset=Post.objects.all().order_by("-date"), template_name='helps/list.html')),
    url(r'^(?P<pk>\d+)$', DetailView.as_view(model=Post, template_name='helps/post.html')),
    url(r'test$', views.viewTest, name='test'),
]