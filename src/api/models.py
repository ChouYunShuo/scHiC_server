from django.db import models
from django.utils import timezone
# Create your models here.


class Dataset(models.Model):
    name = models.CharField(max_length=100)
    file_path = models.TextField(default='')
    create_time = models.DateTimeField(default=timezone.now)
