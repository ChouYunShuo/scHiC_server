from pathlib import Path
import os
import h5py
import numpy as np
import logging
import concurrent.futures

from .utils import parse_region, is_valid_region_string, fetch
logger = logging.getLogger('django')

DATA_DIR = Path(__file__).resolve().parent

FIGURE_WIDTH = 11
FIGURE_HEIGHT = 10


def fetch_spatial(file_path):
    with h5py.File(file_path, 'r') as hdf:
        # Get the res group at the specified root path
        coords = hdf["spatial"].get("coords")
        return np.array(coords, dtype=np.float32)


def fetch_track(file_path, gpath, track_type, res, cell_id, range1):
    with h5py.File(file_path, 'r') as hdf:
        # Get the cell_id group, should have ['pixles', 'indexes', 'ab_score', 'insul_score', 'insulation', 'spatial']
        if track_type not in hdf[gpath]:
            raise ValueError(f"{track_type} is not a valid 1d track type")
        grp = hdf[gpath]
        
        data = grp[f"{track_type}/cell_{cell_id}"][...]
        chrom_offset = grp["indexes/chrom_offset"]
        chrom, start, end = parse_region(range1)

        chr_idx = int(chrom[3:])-1
        start_idx = chrom_offset[chr_idx]+int(start)//res
        end_idx = chrom_offset[chr_idx]+int(end)//res

        temp = [None if np.isnan(x) else float(x)
                for x in data[start_idx:end_idx]]
        return temp


def hic_fetch_map(file_path: str, resolution: str, cell_id: str, range1: str, range2: str):
    # Construct the file path for the hic data and the gpath
    # gpath = os.path.join("resolutions", resolution,
    #                      "cells", f"cell_id{cell_id}")
    gpath = os.path.join("resolutions", resolution)

    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)

    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))

    # Check and get chrom, lo, and hi values for range1 and range2
    if not is_valid_region_string(range1):
        raise ValueError("Request query string not valid: "+str(range1))
    if not is_valid_region_string(range2):
        raise ValueError("Request query string not valid: "+str(range2))

    # Check matrix size

    # Fetch the data from the file path and gpath
    arr = fetch(hic_data_file_path, gpath, cell_id, range1, range2)
    # Apply log2 transformation to the data
    arr = np.log2(arr+1)
    norm_arr = (arr-np.min(arr))/(np.max(arr)-np.min(arr))

    return norm_arr


def fetch_and_add(hic_data_file_path, gpath, id, range1, range2):
    return fetch(hic_data_file_path, gpath, id, range1, range2)


def hic_fetch_maps(file_path: str, resolution: str, cell_ids: list, range1: str, range2: str, isConcurent=False):
    gpath = os.path.join("resolutions", resolution)
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)

    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(hic_data_file_path))

    # Check and get chrom, lo, and hi values for range1 and range2
    if not is_valid_region_string(range1):
        raise ValueError("Request query string not valid: "+str(range1))
    if not is_valid_region_string(range2):
        raise ValueError("Request query string not valid: "+str(range2))

    sum_arr = None

    if (not isConcurent):
        for id in cell_ids:
            arr = fetch(hic_data_file_path, gpath, id, range1, range2)
            if sum_arr is None:
                # If it's the first array, use it to initialize sum_arr
                sum_arr = np.array(arr)
            else:
                # Add the current array to the sum array
                sum_arr += arr
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Use executor.map to parallelize fetch_and_add function
            result_arrays = executor.map(fetch_and_add, [hic_data_file_path]*len(cell_ids), [
                                         gpath]*len(cell_ids), cell_ids, [range1]*len(cell_ids), [range2]*len(cell_ids))

        # Initialize the sum_arr with the first result array
        sum_arr = next(result_arrays)
        for arr in result_arrays:
            sum_arr = np.add(sum_arr, arr)

    sum_arr = np.log2(sum_arr+1)
    norm_arr = (sum_arr-np.min(sum_arr))/(np.max(sum_arr)-np.min(sum_arr))

    return norm_arr


def hic_get_chrom_len(file_path: str, resolution: str, cell_id: str):
    gpath = os.path.join("resolutions", resolution)
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)
    logger.info(hic_data_file_path)
    len_chroms = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(hic_data_file_path))
    with h5py.File(hic_data_file_path, 'r') as hdf:
        # Get the group at the specified root path
        grp = hdf[gpath]
        len_chroms = np.array(grp["chroms"].get("length"))

    return len_chroms


def hic_get_gene_expr(file_path: str, name: str):
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    with h5py.File(hic_data_file_path, 'r') as hdf:
        if name not in hdf["gene_expr"].keys():
            raise IndexError(
                "Gene name {} is invalid.".format(name))
        return np.array(hdf["gene_expr"][name])


def hic_get_meta(file_path: str, meta_type: str):
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)

    meta_vec = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))

    with h5py.File(hic_data_file_path, 'r') as hdf:
        if (meta_type not in list(hdf["meta"])):
            raise ValueError("Request meta type not exist: "+str(meta_type))
        if meta_type == "label":
            vec = hdf["meta"].get("label")
        meta_vec = np.array(vec)
    return meta_vec


def hic_get_embedding(file_path: str, embed_type: str):
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)

    embed_vec = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))

    with h5py.File(hic_data_file_path, 'r') as hdf:
        if embed_type == "pca":
            vec = hdf["embed"].get("pca")
        else:
            vec = hdf["embed"].get("umap")
        embed_vec = np.array(vec, dtype=np.float32)

    return embed_vec


def hic_get_spatial(file_path: str):
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)

    spatial_vec = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    spatial_vec = fetch_spatial(hic_data_file_path)

    return spatial_vec


def hic_get_track(file_path: str, resolution: str, cell_id: str, range1: str, track_type: str):
    # path resolutions/{resolution}/cells/cell_id{cell_id}/track_type
    gpath = os.path.join("resolutions", resolution, "layers/tracks")
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    if not is_valid_region_string(range1):
        raise ValueError("Request query string not valid: "+str(range1))

    return fetch_track(hic_data_file_path, gpath, track_type, int(resolution), cell_id, range1)