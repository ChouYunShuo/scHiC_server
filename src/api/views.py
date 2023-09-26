from . hic.hic_visualizer import hic_fetch_maps, hic_fetch_map, hic_test_api, hic_get_chrom_len, hic_get_embedding, hic_get_spatial, hic_get_track
from django.shortcuts import render, get_object_or_404
from rest_framework.decorators import api_view
from rest_framework.views import APIView
from . models import *
from rest_framework.response import Response
from rest_framework import generics
from .serializers import *

import json
from json import JSONEncoder
import numpy as np
import time
import logging
logger = logging.getLogger('django')


class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.number):
                return obj.tolist()
            elif np.issubdtype(obj.dtype, np.bytes_):
                return obj.astype(np.unicode_).tolist()
        return JSONEncoder.default(self, obj)
# Create your views here.


class scHicQueryView(APIView):
    def post(self, request):
        data = request.data
    
        dataset = get_object_or_404(Dataset, name=data['dataset_name'])
        resolution = data['resolution']
        cell_id = data['cell_id']
        range1 = data['chrom1']
        range2 = data['chrom2']
        logger.info([data['dataset_name'],data['resolution'],data['chrom1'],data['chrom2']])

        try:
            start_time = time.time()
            if(isinstance(cell_id, list)):
                logger.info(len(cell_id))
                arr = hic_fetch_maps(dataset.file_path, resolution, cell_id, range1, range2, True)
            else:
                arr = hic_fetch_map(dataset.file_path, resolution,
                            cell_id, range1, range2)
            end_time = time.time()
            json_array = json.dumps(arr, cls=NumpyArrayEncoder)
            logger.info(f"Runtime: {end_time - start_time} seconds")
            return Response(json_array)
        except Exception as e:
            return Response({"invalid": str(e)}, status=400)


@api_view(["POST"])
def chromLenView(request):
    data = request.data
    logger.info(data)

    resolution = data['resolution']
    cell_id = data['cell_id']
    dataset = get_object_or_404(Dataset, name=data['name'])
    try:
        arr = hic_get_chrom_len(dataset.file_path, resolution, cell_id)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        #logger.info(json_array)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)
    
@api_view(["POST"])
def embeddingView(request):
    data = request.data
    logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        resolution = data['resolution']
        embed_type = data['embed_type']
        arr = hic_get_embedding(dataset.file_path, resolution, embed_type)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)

@api_view(["POST"])
def trackView(request):
    data = request.data
    logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        resolution = data['resolution']
        track_type = data['type']
        cell_id = data['cell_id']
        range1 = data['chrom1']
        if(isinstance(cell_id, list)):
            arr = hic_get_track(dataset.file_path, resolution, cell_id[0], range1, track_type)
        else:
            arr = hic_get_track(dataset.file_path, resolution, cell_id, range1, track_type)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)

@api_view(["POST"])
def spatialView(request):
    data = request.data
    logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        resolution = data['resolution']
        gene_name = data['gene_name']
        arr = hic_get_spatial(dataset.file_path, resolution, gene_name)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)

class scHicTestView(APIView):
    def post(self, request):
        data = request.data
        logger.info(data)

        resolution = data['resolution']
        range1 = data['chrom1']
        range2 = data['chrom2']
        try:
            arr = hic_test_api(resolution, range1, range2)
            json_array = json.dumps(arr, cls=NumpyArrayEncoder)
            return Response(json_array)
        except Exception as e:
            return Response({"invalid": str(e)}, status=400)


class DatasetListAPIView(generics.ListAPIView):
    queryset = Dataset.objects.all()
    serializer_class = DatasetSerializer


@api_view(["GET"])
def dataset_retreive_view(request, pk=None, *args, **kwargs):
    if pk is not None:
        obj = get_object_or_404(Dataset, pk=pk)
        data = DatasetSerializer(obj, many=False).data
        return Response(data)
    queryset = Dataset.objects.all()
    data = DatasetSerializer(queryset, many=True).data
    return Response(data)


"""
{
"file_path": "data/scHiC.h5",
"grp_path": "resolutions/100000/cells/cell_id0",
"range1": "chrom2:0-40000000",
"range2": "chrom2:0-40000000"
}
"""
