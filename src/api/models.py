from django.db import models
from django.utils import timezone
# Create your models here.


class Dataset(models.Model):
    name = models.CharField(max_length=100)
    file_path = models.TextField(default='')
    description = models.TextField(blank=True, null=True)
    resolutions = models.TextField(default='500000')
    cells = models.IntegerField(default=4)
    create_time = models.DateTimeField(default=timezone.now)
