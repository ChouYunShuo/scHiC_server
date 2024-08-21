from . hic.hic_visualizer import hic_fetch_maps,hic_fetch_group, hic_fetch_map, hic_get_chrom_len, hic_get_embedding, hic_get_spatial, hic_get_track, hic_get_meta, hic_get_gene_expr, hic_get_session_json,hic_upload_session_json
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
        cell_id = data.get('cell_id', None)
        cell_type = data.get('cell_type', None)  
        range1 = data['chrom1']
        range2 = data['chrom2']
        logger.info([data['dataset_name'], resolution,cell_id,cell_type,
                    range1, range2])

        try:
            start_time = time.time()
            if cell_type != "":
                arr = hic_fetch_group(
                    dataset.file_path, resolution, cell_type, range1, range2)
            elif (isinstance(cell_id, list)):
                logger.info(len(cell_id))
                arr = hic_fetch_maps(
                    dataset.file_path, resolution, cell_id, range1, range2, True)
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
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)


@api_view(["POST"])
def embeddingView(request):
    data = request.data
    # logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        embed_type = data['embed_type']
        arr = hic_get_embedding(dataset.file_path, embed_type)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)


@api_view(["POST"])
def metaView(request):
    data = request.data
    # logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        meta_type = data['meta_type']
        arr = hic_get_meta(dataset.file_path, meta_type)
        logger.info(arr)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)


@api_view(["POST"])
def spatialView(request):
    data = request.data
    # logger.info(data)
    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        arr = hic_get_spatial(dataset.file_path)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)


@api_view(["POST"])
def geneExprView(request):
    data = request.data
    # logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        name = data['name']
        arr = hic_get_gene_expr(dataset.file_path, name)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)


@api_view(["POST"])
def trackView(request):
    data = request.data
    # logger.info(data)

    dataset = get_object_or_404(Dataset, name=data['dataset_name'])
    try:
        resolution = data['resolution']
        track_type = data['type']
        cell_id = data['cell_id']
        range1 = data['chrom1']
        if (isinstance(cell_id, list)):
            arr = hic_get_track(dataset.file_path, resolution,
                                cell_id[0], range1, track_type)
        else:
            arr = hic_get_track(dataset.file_path, resolution,
                                cell_id, range1, track_type)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)

@api_view(["GET"])
def session_retreive_view(request, uuid):
    obj = get_object_or_404(Session, pk=uuid)
    data = SessionSerializer(obj, many=False).data
    file_path = data['file_path']
    logger.info(uuid,file_path)
    session_json = hic_get_session_json(file_path)
    return Response(session_json)

@api_view(["POST"])
def session_upload_view(request):
    data = request.data
    try:
        config = data['config']
        file_name = data['file_name']
        logger.info(config, file_name)
        session_uuid = hic_upload_session_json(config, file_name)
        return Response({"session_uuid": session_uuid})
    except Exception as e:
        return Response({"invalid": str(e)}, status=400)

    


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
