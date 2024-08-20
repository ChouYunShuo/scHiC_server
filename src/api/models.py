from django.db import models
from django.utils import timezone
import uuid
# Create your models here.


class Dataset(models.Model):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, unique=True)
    session_uuid = models.UUIDField(null=True, blank=True)
    name = models.CharField(max_length=100)
    file_path = models.TextField(default='')
    description = models.TextField(blank=True, null=True)
    resolutions = models.TextField(default='500000')
    cells = models.IntegerField(default=4)
    create_time = models.DateTimeField(default=timezone.now)

class Session(models.Model):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, unique=True)
    name = models.CharField(max_length=100)
    file_path = models.TextField(default='')
    create_time = models.DateTimeField(default=timezone.now)