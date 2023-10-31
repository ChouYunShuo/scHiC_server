from pathlib import Path
import os
from typing import Tuple
from django.conf import settings
import base64
from io import BytesIO
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import h5py
from scipy.sparse import coo_matrix
import numpy as np
import matplotlib
import re
import logging
import hicstraw
import time
import concurrent.futures
from .utils import rlencode #need . to import from same dir

logger = logging.getLogger('django')
matplotlib.use('Agg')

#ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
DATA_DIR = Path(__file__).resolve().parent
REDMAP = LinearSegmentedColormap.from_list(
    "bright_red", [(1, 1, 1), (1, 0, 0)])


def parse_region(range: str) -> Tuple:
    """
        range: chrom1:10000000-20000000
    """

    # future1: do checking on valid string

    parts = range.split(":")
    pos = parts[1].split("-")

    return parts[0], pos[0], pos[1]


def is_valid_region_string(region_string):
    """
    Check if a string is in the format "chromX:start-end" where X is a number between 1 and 23 or X,Y and start < end, end is greater than 0
    :param region_string: The input string to check
    :return: A boolean indicating whether the string is in the correct format
    """
    pattern = re.compile(
        r"^chrom([1-9]|1[0-9]|2[0-3]|X|Y):([0-9][0-9]*)-([1-9][0-9]*)$")
    match = pattern.match(region_string)
    if match:
        start, end = float(match.group(2)), float(match.group(3))
        return start < end
    return False


def get_chrom_ids(name_chroms):

    chrom_ids = {}
    id = 0
    for n in name_chroms:
        chrom_ids["chrom"+n.astype(str)] = id
        id += 1

    return chrom_ids


def get_chrom_dict(name_chroms, size_chroms):

    chrom_dict = {}
    for n, s in zip(name_chroms, size_chroms):
        chrom_dict["chrom"+n.astype(str)] = s

    return chrom_dict


def region_to_extent(res_grp, h5, chrom_ids, region):
    """
    Convert a region to the corresponding extent in an HDF5 file.
    :param h5: The HDF5 file object.
    :param chrom_ids: A dictionary containing the chromosome IDs.
    :param region: A tuple of the form (chromosome, start, end)
    :return: A tuple containing the starting and ending indices of the extent.
    """
    chrom, start, end = region
    cid = chrom_ids[chrom]

    chrom_lo = h5["indexes"]["chrom_offset"][cid]
    chrom_hi = h5["indexes"]["chrom_offset"][cid + 1]
    chrom_bins = res_grp["bins"]["start"][chrom_lo:chrom_hi]

    return chrom_lo + np.searchsorted(chrom_bins, int(start), "right") - 1, chrom_lo + np.searchsorted(chrom_bins, int(end), "left")

def fetch(file_path, root_path, cell_id, region1, region2):
    """
    Retrieve data from an HDF5 file for a specified pair of regions.
    :param file_path: The path to the HDF5 file.
    :param root_path: The root path of the group containing the data.
    :param region1: The first region to retrieve data for.
    :param region2: The second region to retrieve data for.
    :return: A dense array containing the data for the specified regions.
    """
     # Open an HDF5 file for reading using h5py
    with h5py.File(file_path, 'r') as hdf:
        # Get the group at the specified root path
        grp = hdf[root_path]
        cell_grp = grp.get("cells/cell_id"+str(cell_id))
        # Get the array of chromosome names
        name_chroms = np.array(grp["chroms"].get("name"))
        # Get the chromosome IDs for each region
        chromids = get_chrom_ids(name_chroms)
        
        region1 = parse_region(region1)
        region2 = parse_region(region2)

        # Get the extent of the first region in the group
        i0, i1 = region_to_extent(grp, cell_grp, chromids, region1)
        # Get the extent of the second region in the group
        j0, j1 = region_to_extent(grp, cell_grp, chromids, region2)

        # Perform a query on the group to get the values at the specified extent
        i, j, v = query(cell_grp, i0, i1, j0, j1)

        # Check if i, j, and v have data
        if not i.size or not j.size or not v.size:
            raise ValueError("No data found for the specified regions.")

        # Return a dense array created from the sparse coo_matrix
        arr = coo_matrix((v, (i - i0, j - j0)), shape=(i1 - i0, j1 - j0)).toarray()
        return arr



def query(h5, i0, i1, j0, j1, field="count") -> Tuple:
    """
    Retrieve data from an HDF5 file for a specified range of rows and columns.
    :param h5: The HDF5 file object.
    :param i0: The starting row index.
    :param i1: The ending row index.
    :param j0: The starting column index.
    :param j1: The ending column index.
    :param field: The field of the data to retrieve (default is "count").
    :return: A tuple of arrays containing the row indices, column indices, and values of the data.
    """
    i, j, v = [], [], []

    edges = h5["indexes"]["bin1_offset"][i0: i1 + 1]
    data = h5["pixels"][field]
    p0, p1 = edges[0], edges[-1]
    all_bin2 = h5["pixels"]["bin2_id"][p0:p1]
    all_data = data[p0:p1]
    dtype = all_bin2.dtype
    for row_id, lo, hi in zip(
        range(i0, i1), edges[:-1] - p0, edges[1:] - p0
    ):
        bin2 = all_bin2[lo:hi]
        mask = (bin2 >= j0) & (bin2 < j1)
        cols = bin2[mask]

        i.append(np.full(len(cols), row_id, dtype=dtype))
        j.append(cols)
        v.append(all_data[lo:hi][mask])

    i = np.concatenate(i, axis=0)
    j = np.concatenate(j, axis=0)
    v = np.concatenate(v, axis=0)

    return i, j, v


FIGURE_WIDTH = 11
FIGURE_HEIGHT = 10

def fetch_spatial(file_path,gpath,gene_name):
    with h5py.File(file_path, 'r') as hdf:
        # Get the res group at the specified root path
        grp = hdf[gpath]
        
        coords = grp["spatial"].get("coords") 
        gene_expr =  grp["spatial"].get(gene_name)
        spatial_vec = np.column_stack((np.array(coords, dtype=np.float32), np.array(gene_expr)))

        return spatial_vec
def fetch_insulation(file_path, gpath, res, cell_id, range1):
    with h5py.File(file_path, 'r') as hdf:
        # Get the resolution group, should have ['bins', 'cells', 'chroms', 'embed', 'insulation', 'spatial']
        grp = hdf[gpath]
        if "insulation" not in grp:
            raise ValueError("track type: insulation does not exist in dataset")
        
        n_chroms = len(np.unique(grp["insulation/bin/chrom"]))
        n_bins = len(grp["insulation/bin/chrom"])
        chrom_ids = grp["insulation/bin/chrom"]
        data = grp["insulation/cell_"+str(cell_id)][...]
        chrom_offset = np.zeros(n_chroms + 1, dtype=np.int32)
        index = 0
        for start, length, value in zip(*rlencode(chrom_ids)):
            chrom_offset[index] = start
            index += 1
        chrom_offset[index] = n_bins

        chrom, start ,end = parse_region(range1)
        
        chr_idx = int(chrom[5:])-1
        start_idx = chrom_offset[chr_idx]+int(start)//res
        end_idx = chrom_offset[chr_idx]+int(end)//res
    
        temp = [None if np.isnan(x) else x for x in data[start_idx:end_idx]]
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
    
    if(not isConcurent):
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
            result_arrays = executor.map(fetch_and_add, [hic_data_file_path]*len(cell_ids), [gpath]*len(cell_ids), cell_ids, [range1]*len(cell_ids), [range2]*len(cell_ids))

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
   
    len_chroms = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(hic_data_file_path))
    with h5py.File(hic_data_file_path, 'r') as hdf:
        # Get the group at the specified root path
        grp = hdf[gpath]
        len_chroms = np.array(grp["chroms"].get("length"))

    return len_chroms

def hic_get_embedding(file_path: str, resolution: str, embed_type: str):
    gpath = os.path.join("resolutions", resolution)
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)
    
   
    embed_vec = []
    if not os.path.exists(hic_data_file_path):
        
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    
    with h5py.File(hic_data_file_path, 'r') as hdf:
        
        # Get the res group at the specified root path
        grp = hdf[gpath]
        if embed_type== "pca":
            vec = grp["embed"].get("pca") 
        else:
            vec = grp["embed"].get("umap") 

        cell_type =  grp["embed"].get("label")

        embed_vec = np.column_stack((np.array(vec, dtype=np.float32), np.array(cell_type)))
    return embed_vec
        
def hic_get_spatial(file_path: str, resolution: str, gene_name: str):
    gpath = os.path.join("resolutions", resolution)
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)
   
    spatial_vec = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    spatial_vec = fetch_spatial(hic_data_file_path, gpath,gene_name)

    return spatial_vec

def hic_get_track(file_path: str, resolution: str, cell_id: str, range1: str, track_type: str):
    gpath = os.path.join("resolutions", resolution)
    hic_data_file_path = os.path.join(DATA_DIR, "hic_data", file_path)
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    if not is_valid_region_string(range1):
        raise ValueError("Request query string not valid: "+str(range1))
    if track_type == 'insulation':
        return fetch_insulation(hic_data_file_path, gpath, int(resolution), cell_id,range1)
    
    return []

def hic_test_api(resolution: str, range1: str, range2: str):
    hic_path = os.path.join(DATA_DIR, "hic_data", "ENCFF718AWL.hic")
    hic = hicstraw.HiCFile(hic_path)
    if not is_valid_region_string(range1):
        raise ValueError("Request query string not valid: "+str(range1))
    if not is_valid_region_string(range2):
        raise ValueError("Request query string not valid: "+str(range2))
    row_chrom, row_lo, row_hi = parse_region(range1)
    col_chrom, col_lo, col_hi = parse_region(range2)

    mzd = hic.getMatrixZoomData(
        row_chrom[5:], col_chrom[5:], "observed", "KR", "BP", int(resolution))
    arr = mzd.getRecordsAsMatrix(
        int(row_lo), int(row_hi), int(col_lo), int(col_hi))
    arr = np.log2(arr+1)
    arr = (arr-np.min(arr))/(np.max(arr)-np.min(arr))

    return arr


def visualize(fpath, range1, range2):
    #c = Hic_wrapper(cool_uri)

    # future 1: check if use chroms size for end

    row_chrom, row_lo, row_hi = parse_region(range1)
    col_chrom, col_lo, col_hi = parse_region(range2)

    # future 2: check matrix size
    #file_path, root_path, cell_id
    cnt = 1
    cell_ids = [str(i) for i in range(cnt)]

    start_time = time.time()
    arr = hic_fetch_maps(fpath,"100000",cell_ids,range1, range2, True)
    end_time = time.time()
    print(f"Runtime: {end_time - start_time} seconds")

    plt.figure(figsize=(11, 10))
    plt.title(range1 + "&" + range2)
    plt.imshow(
        arr,
        interpolation="none",
        extent=[int(col_lo), int(col_hi), int(row_hi), int(row_lo)],
        cmap=REDMAP,
    )
    plt.savefig(f"mygraph_{cnt}.png")


if __name__ == "__main__":
    
    # fpath = os.path.join(ROOT_DIR, "hic_data", "4DN_scHi-C_Kim.h5")
    # gpath = os.path.join("resolutions", "500000")
    # visualize(fpath, gpath, "chrom1:0-243199372", "chrom1:0-243199372")
    
    #hic_get_chrom_len('scHiC5.h5', '50000', '0')

    root_path = Path(__file__).resolve().parent.parent.parent.parent
    file_path = str(root_path)+"/hic_data/mouse2_slice99.h5"
    
    visualize(file_path,"chrom1:0-19515427", "chrom1:0-19515427")
    """
    Test for Embedding
    with h5py.File(file_path, 'r') as hdf:
            all_res = hdf.get("resolutions")
            print(list(all_res))
            cur_res = hdf.get("resolutions/500000")

            vec_pca = cur_res["embed"].get("pca") 
            vec_umap = cur_res["embed"].get("umap") 
            cell_type =  cur_res["embed"].get("label")

            fig = plt.figure(figsize=(14, 5))
            ax = plt.subplot(1, 2, 1)
            sns.scatterplot(x=vec_pca[:, 0], y=vec_pca[:, 1], hue=cell_type, ax=ax, s=6, linewidth=0)
            handles, labels = ax.get_legend_handles_labels()
            labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
            ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
            ax = plt.subplot(1, 2, 2)
        
            sns.scatterplot(x=vec_umap[:, 0], y=vec_umap[:, 1], hue=cell_type, ax=ax, s=6, linewidth=0)
            handles, labels = ax.get_legend_handles_labels()
            labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
            ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
            plt.tight_layout()
            plt.savefig("embed.png")
    """
    