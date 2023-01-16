from . hic.hic_visualizer import hic_api
from django.shortcuts import render, get_object_or_404
from rest_framework.views import APIView
from . models import *
from rest_framework.response import Response
from rest_framework import generics
from .serializers import *

import json
from json import JSONEncoder
import numpy
import logging
logger = logging.getLogger('django')


class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)
# Create your views here.


class scHicQueryView(APIView):
    def post(self, request):
        data = request.data
        logger.info(data)

        dataset = get_object_or_404(Dataset, name=data['dataset_name'])
        resolution = data['resolution']
        cell_id = data['cell_id']
        range1 = data['chrom1']
        range2 = data['chrom2']
        try:
            arr = hic_api(dataset.file_path, resolution,
                          cell_id, range1, range2)
            json_array = json.dumps(arr, cls=NumpyArrayEncoder)
            return Response(json_array)
        except Exception as e:
            return Response({"invalid": str(e)}, status=400)


class DatasetListAPIView(generics.ListAPIView):
    queryset = Dataset.objects.all()
    serializer_class = DatasetSerializer


"""
{
"file_path": "data/scHiC.h5",
"grp_path": "resolutions/100000/cells/cell_id0",
"range1": "chrom2:0-40000000",
"range2": "chrom2:0-40000000"
}
"""
