from typing import Tuple
import h5py
from scipy.sparse import coo_matrix
import numpy as np
import re

def parse_region(range: str) -> Tuple:
    """
        range: chr1:10000000-20000000
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
        r"^chr([1-9]|1[0-9]|2[0-3]|X|Y):([0-9][0-9]*)-([1-9][0-9]*)$")
    match = pattern.match(region_string)
    if match:
        start, end = float(match.group(2)), float(match.group(3))
        return start < end
    return False


def get_chrom_ids(name_chroms):

    chrom_ids = {}
    id = 0
    for n in name_chroms:
        chrom_ids["chr"+n.astype(str)] = id
        id += 1

    return chrom_ids


def get_chrom_dict(name_chroms, size_chroms):

    chrom_dict = {}
    for n, s in zip(name_chroms, size_chroms):
        chrom_dict["chr"+n.astype(str)] = s

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
        # FIXME
        cell_grp = grp.get("layers/imputed_0neighbor/cell_"+str(cell_id))
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
        arr = coo_matrix((v, (i - i0, j - j0)),
                         shape=(i1 - i0, j1 - j0)).toarray()
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

    # edges = h5["indexes"]["bin1_offset"][i0: i1 + 1]
    edges = h5["indexes"]["bin1_offset"][i0: i1 + 1]
    data = h5["pixels"][field]  # pixels
    p0, p1 = edges[0], edges[-1]
    all_bin2 = h5["pixels"]["bin2_id"][p0:p1]  # pixels
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