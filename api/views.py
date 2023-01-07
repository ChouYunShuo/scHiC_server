from . hic.hic_visualizer import hic_api
from django.shortcuts import render
from rest_framework.views import APIView
from . models import *
from rest_framework.response import Response
from . hic_serializer import *

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


class scHicView(APIView):
    def post(self, request):
        data = request.data
        logger.info(data)
        range1 = data['chrom1']
        range2 = data['chrom2']
        arr = hic_api(range1, range2)
        json_array = json.dumps(arr, cls=NumpyArrayEncoder)
        return Response(json_array)


"""
{
"file_path": "data/scHiC.h5",
"grp_path": "resolutions/100000/cells/cell_id0",
"range1": "chrom2:0-40000000",
"range2": "chrom2:0-40000000"
}
"""
