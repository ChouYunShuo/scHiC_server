from django.contrib import admin
from .models import Dataset, Session
# Register your models here.

@admin.register(Dataset)
class DatasetAdmin(admin.ModelAdmin):
    list_display = ('uuid', 'session_uuid', 'name', 
                    'file_path','description', 'resolutions', 'cells'
    )

@admin.register(Session)
class SessionAdmin(admin.ModelAdmin):
    list_display = ('uuid', 'name', 'file_path',)  # Adjust the fields to display
