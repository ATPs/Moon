from django.db import models

# Create your models here.
class Post(models.Model):
    title = models.CharField(max_length=140)
    body = models.TextField()
    date = models.DateTimeField()
    author = models.CharField(max_length=140,default="XC")