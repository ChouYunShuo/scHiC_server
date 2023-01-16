from pathlib import Path
from django.conf import settings
import base64
from io import BytesIO
import os
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import h5py
from scipy.sparse import coo_matrix
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
#from . import Hic_wrapper

ROOT_DIR = Path(__file__).resolve().parent.parent.parent.parent
REDMAP = LinearSegmentedColormap.from_list(
    "bright_red", [(1, 1, 1), (1, 0, 0)])


def parse_region(range: str, chrom_dict=None):
    """
        range: chrom1:10000000-20000000
    """

    # future1: do checking on valid string

    parts = range.split(":")
    pos = parts[1].split("-")

    return (parts[0], pos[0], pos[1])


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
    def _region_to_extent():
        chrom, start, end = region
        cid = chrom_ids[chrom]

        chrom_lo = h5["indexes"]["chrom_offset"][cid]
        chrom_hi = h5["indexes"]["chrom_offset"][cid + 1]
        chrom_bins = h5["bins"]["start"][chrom_lo:chrom_hi]

        yield chrom_lo + np.searchsorted(chrom_bins, int(start), "right") - 1
        yield chrom_lo + np.searchsorted(chrom_bins, int(end), "left")

    return tuple(_region_to_extent())


def fetch(file_path, root_path, region1, region2):
    with h5py.File(file_path, 'r') as hdf:
        grp = hdf[root_path]

        name_chroms = np.array(grp["chroms"].get("name"))
        #chrom_dict = get_chrom_dict(name_chroms, size_chroms)
        # print(chrom_dict)

        region1 = parse_region(region1)
        region2 = parse_region(region2)

        chromids = get_chrom_ids(name_chroms)

        i0, i1 = region_to_extent(grp, chromids, region1)
        j0, j1 = region_to_extent(grp, chromids, region2)

        i, j, v = query(grp, i0, i1, j0, j1)

        return coo_matrix((v, (i - i0, j - j0)), shape=(i1 - i0, j1 - j0)).toarray()


def query(h5, i0, i1, j0, j1, field="count"):
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


def hic_api(fp, resolution, cell_id, range1, range2):
    gpath = "resolutions/"+resolution+"/cells/cell_id"+cell_id
    fpath = str(ROOT_DIR)+"/hic_data/"+fp

    row_chrom, row_lo, row_hi = parse_region(range1)
    col_chrom, col_lo, col_hi = parse_region(range2)

    # future 2: check matrix size

    arr = fetch(fpath, gpath, range1, range2)
    arr = np.log2(arr+1)

    image = BytesIO()
    im = plt.figure(figsize=(11, 10))
    plt.title(range1 + "&" + range2)
    plt.imshow(
        arr,
        interpolation="none",
        extent=[int(col_lo), int(col_hi), int(row_hi), int(row_lo)],
        cmap=REDMAP,
    )
    im.savefig(image, format='png')
    image.seek(0)  # rewind to beginning of file
    return base64.b64encode(image.getvalue()).decode('utf-8')


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
    pass
