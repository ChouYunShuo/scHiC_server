from rest_framework import serializers
from .models import Dataset, Session


class DatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dataset
        fields = [
            'uuid',
            'session_uuid',
            'name',
            'description',
            'resolutions',
            'cells'
        ]
class SessionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Session
        fields = [
            'uuid',
            'name',
            'file_path',
        ]
