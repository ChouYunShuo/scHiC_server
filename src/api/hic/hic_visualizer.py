from pathlib import Path
import os
from typing import Tuple
from django.conf import settings
import base64
from io import BytesIO
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import h5py
from scipy.sparse import coo_matrix
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
import logging
import hicstraw

logger = logging.getLogger('django')
matplotlib.use('Agg')

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
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


def region_to_extent(h5, chrom_ids, region):
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
    chrom_bins = h5["bins"]["start"][chrom_lo:chrom_hi]

    return chrom_lo + np.searchsorted(chrom_bins, int(start), "right") - 1, chrom_lo + np.searchsorted(chrom_bins, int(end), "left")


def fetch(file_path, root_path, region1, region2):
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

        # Get the array of chromosome names
        name_chroms = np.array(grp["chroms"].get("name"))
        # Get the array of chromosome length
        len_chroms = np.array(grp["chroms"].get("length"))
        # Get the chromosome IDs for each region
        chromids = get_chrom_ids(name_chroms)

        region1 = parse_region(region1)
        region2 = parse_region(region2)

        # Get the extent of the first region in the group
        i0, i1 = region_to_extent(grp, chromids, region1)
        # Get the extent of the second region in the group
        j0, j1 = region_to_extent(grp, chromids, region2)

        # Perform a query on the group to get the values at the specified extent
        i, j, v = query(grp, i0, i1, j0, j1)

        # Check if i, j, and v have data
        if not i.size or not j.size or not v.size:
            raise ValueError("No data found for the specified regions.")

        # Return a dense array created from the sparse coo_matrix
        return coo_matrix((v, (i - i0, j - j0)), shape=(i1 - i0, j1 - j0)).toarray()


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


def hic_api(file_path: str, resolution: str, cell_id: str, range1: str, range2: str):
    # Construct the file path for the hic data and the gpath
    gpath = os.path.join("resolutions", resolution,
                         "cells", f"cell_id{cell_id}")
    hic_data_file_path = os.path.join(ROOT_DIR, "hic_data", file_path)
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))

    # Check and get chrom, lo, and hi values for range1 and range2
    if not is_valid_region_string(range1):
        raise ValueError("Request query string not valid: "+str(range1))
    if not is_valid_region_string(range2):
        raise ValueError("Request query string not valid: "+str(range2))
    row_chrom, row_lo, row_hi = parse_region(range1)
    col_chrom, col_lo, col_hi = parse_region(range2)

    # Check matrix size

    # Fetch the data from the file path and gpath
    arr = fetch(hic_data_file_path, gpath, range1, range2)
    # Apply log2 transformation to the data
    arr = np.log2(arr+1)
    norm_arr = (arr-np.min(arr))/(np.max(arr)-np.min(arr))

    return norm_arr


def hic_get_chrom_len(file_path: str, resolution: str, cell_id: str):
    gpath = os.path.join("resolutions", resolution,
                         "cells", f"cell_id{cell_id}")
    hic_data_file_path = os.path.join(ROOT_DIR, "hic_data", file_path)
    len_chroms = []
    if not os.path.exists(hic_data_file_path):
        raise FileNotFoundError(
            "No such file or directory with name: "+str(file_path))
    with h5py.File(hic_data_file_path, 'r') as hdf:
        # Get the group at the specified root path
        grp = hdf[gpath]
        len_chroms = np.array(grp["chroms"].get("length"))

    return len_chroms


def hic_test_api(resolution: str, range1: str, range2: str):
    root_dir = Path(__file__).resolve().parent.parent.parent.parent
    hic_path = str(root_dir)+"/hic_data/ENCFF718AWL.hic"
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


def visualize(fpath, gpath, range1, range2):
    #c = Hic_wrapper(cool_uri)

    # future 1: check if use chroms size for end

    row_chrom, row_lo, row_hi = parse_region(range1)
    col_chrom, col_lo, col_hi = parse_region(range2)

    # future 2: check matrix size

    arr = fetch(fpath, gpath, range1, range2)
    arr = np.log2(arr+1)
    plt.figure(figsize=(11, 10))
    plt.title(range1 + "&" + range2)
    plt.imshow(
        arr,
        interpolation="none",
        extent=[int(col_lo), int(col_hi), int(row_hi), int(row_lo)],
        cmap=REDMAP,
    )
    plt.savefig("mygraph.png")


if __name__ == "__main__":
    '''
    fpath = os.path.join(ROOT_DIR, "hic_data", "scHiC5.h5")
    gpath = os.path.join("resolutions", "50000",
                         "cells", f"cell_id{0}")
    visualize(fpath, gpath, "chrom1:0-20000000", "chrom1:0-20000000")
    '''
    hic_get_chrom_len('scHiC5.h5', '50000', '0')
