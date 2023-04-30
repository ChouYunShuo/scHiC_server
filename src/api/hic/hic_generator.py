import time
import numpy as np
import h5py
from collections import defaultdict
from utils import rlencode
import seaborn as sns
import hicstraw
import os
import logging
import traceback
from pathlib import Path

CHROM_DTYPE = np.dtype("S")
CHROMID_DTYPE = np.int32
CHROMSIZE_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32
OFFSET_DTYPE = np.int64

"""
└── resolutions
     ├── 1000
     │   ├── bins
     │   ├── chroms
     │   └── cells
     │       ├── cell_id1
     │       │   ├── bins
     │       │   ├── chroms
     │       │   ├── pixels
     │       │   └── indexes
     │       ├── cell_id2
     │       │   ├── bins
     │       │   ├── chroms
     │       │   ├── pixels
     │       │   └── indexes
     ├── 5000
     │   ├── bins
     │   ├── chroms
     │   └── cells
     │       ├── cell_id1
     │       │   ├── bins
     │       │   ├── chroms
     │       │   ├── pixels
     │       │   └── indexes
     │       ├── cell_id2
     │       │   ├── bins
     │       │   ├── chroms
     │       │   ├── pixels
     │       │   └── indexes
     ├── 5000
     │   
     ├── 10000
     │   
"""


def check_resolution_valid(res, hic_res):
    if not res:
        res = hic_res
    else:
        check = all(item in hic_res for item in res)
        if not check:
            raise RuntimeError(
                "Provided resolutions not provided in .hic file")

    return res


def write_chroms(grp, hic, h5_opts):
    chrom_dict = defaultdict(list)
    for chrom in hic.getChromosomes():
        if chrom.name == 'All':
            continue
        chrom_dict["name"].append(chrom.name)
        chrom_dict["length"].append(chrom.length)

    n_chroms = len(chrom_dict["name"])

    names = np.array(chrom_dict["name"], dtype=CHROM_DTYPE)
    lengths = np.array(chrom_dict["length"], dtype=CHROMSIZE_DTYPE)

    grp.create_dataset('name', shape=(n_chroms,),
                       dtype=names.dtype, data=names, **h5_opts)
    grp.create_dataset("length", shape=(n_chroms,),
                       dtype=lengths.dtype, data=lengths, **h5_opts)

    return names, lengths


def write_bins(res, chroms_names, chrom_lens, grp, h5_opts):
    bin_dict = defaultdict(list)
    for index in range(len(chroms_names)):
        for i in range(int((chrom_lens[index])/(res))+1):
            bin_dict["chrom"].append(chroms_names[index])
            bin_dict["start"].append(i*res)
            bin_dict["end"].append((i+1)*res)

    n_bins = len(bin_dict["chrom"])

    chroms = np.array(bin_dict["chrom"], dtype=CHROM_DTYPE)
    starts = np.array(bin_dict["start"], dtype=COORD_DTYPE)
    ends = np.array(bin_dict["end"], dtype=COORD_DTYPE)

    grp.create_dataset('chrom', shape=(n_bins,),
                       dtype=chroms.dtype, data=chroms, **h5_opts)
    grp.create_dataset("start", shape=(n_bins,),
                       dtype=starts.dtype, data=starts, **h5_opts)
    grp.create_dataset("end", shape=(n_bins,),
                       dtype=ends.dtype, data=ends, **h5_opts)


def setup_pixels(grp, nbins, h5_opts):
    max_shape = nbins*nbins
    grp.create_dataset('bin1_id', shape=(max_shape,),
                       dtype=BIN_DTYPE, **h5_opts)
    grp.create_dataset("bin2_id", shape=(max_shape,),
                       dtype=BIN_DTYPE, **h5_opts)
    grp.create_dataset("count", shape=(max_shape,),
                       dtype=COUNT_DTYPE, **h5_opts)


def write_pixels(grp, hic, chrom_names, chrom_lens, chunksize, res, columns, chrom_offset):
    cellMatrixParser = MatrixParser(
        hic, chrom_names, chrom_lens, chunksize, res, chrom_offset)
    m_size = 0
    for chunk in cellMatrixParser:
        dsets = [grp[col] for col in columns]
        n = len(chunk[columns[0]])
        for col, dset in zip(columns, dsets):
            dset.resize((m_size + n,))
            dset[m_size: m_size + n] = chunk[col]
        m_size += n


def sort_pixel_by_bin(grp, cols):
    dsets = [np.array(grp[col]) for col in cols]
    ind = np.lexsort((dsets[1], dsets[0]))

    for i, col in enumerate(cols):
        grp[col][...] = dsets[i][ind]

# with h5py.File("data/scHiC2.h5", 'r+') as hdf:
   ## G1 = hdf.get("resolutions/100000/cells/cell_id0")
   # grp = G1.get("pixels")
   # sort_pixel_by_bin(grp, list(grp))


def get_bin_index(grp, n_chroms, n_bins):
    chrom_ids = grp["chrom"]
    chrom_offset = np.zeros(n_chroms + 1, dtype=OFFSET_DTYPE)
    index = 0
    for start, length, value in zip(*rlencode(chrom_ids)):
        chrom_offset[index] = start
        index += 1
    chrom_offset[index] = n_bins

    return chrom_offset


def get_pixel_index(grp, n_bins, n_pixels):
    bin1 = np.array(grp["bin1_id"])
    bin1_offset = np.zeros(n_bins + 1, dtype=OFFSET_DTYPE)
    curr_val = 0

    for start, length, value in zip(*rlencode(bin1, 1000000)):
        bin1_offset[curr_val: value + 1] = start
        curr_val = value+1

    bin1_offset[curr_val:] = n_pixels

    return bin1_offset


def write_index(grp, chrom_offset, bin_offset, h5_opts):

    grp.create_dataset(
        "chrom_offset",
        shape=(len(chrom_offset),),
        dtype=OFFSET_DTYPE,
        data=chrom_offset, **h5_opts
    )
    grp.create_dataset(
        "bin1_offset",
        shape=(len(bin_offset),),
        dtype=OFFSET_DTYPE,
        data=bin_offset, **h5_opts
    )


'''

Helper class that transform .hic to hdf5 format.

'''


class SCHiCGenerator:
    def __init__(self, hic_path, multi_processed=True):
        self.hic_path = hic_path
        self.hic = hicstraw.HiCFile(hic_path)
        self.h5_path = ""
        self.h5_opts = {}

    def set_h5_opts(self):
        self.h5_opts["compression"] = "gzip"
        self.h5_opts["compression_opts"] = 6
        self.h5_opts["shuffle"] = True
        self.h5_opts["chunks"] = True

    def create_all_h5(self, file_name, res=[]):

        if os.path.exists(file_name):
            raise RuntimeError("sc-HiC file: " + file_name + " already exists")
        self.h5_path = file_name
        self.set_h5_opts()
        all_res = check_resolution_valid(res, self.hic.getResolutions())
        with h5py.File(self.h5_path, 'w') as hdf:
            root_grp = hdf.create_group("resolutions")
        for res in all_res:
            print("Creating resolution: "+str(res))
            self.create_res_h5(root_grp, res, self.h5_opts)

    '''
    params:
        res: current resolution
        cell_cnt: numbers of single cell
        chunksize: bin size for query, endx = startx+res*chunksize
    '''

    def create_res_h5(self, root_grp, res, h5_opts, cell_cnt=4, chunksize=500):

        with h5py.File(self.h5_path, 'r+') as hdf:
            res_grp = hdf.create_group("resolutions/"+str(res))

            grp_chroms = res_grp.create_group("chroms")
            np_chroms_names, np_chroms_length = write_chroms(
                grp_chroms, self.hic, h5_opts)

            grp_bins = grp_chroms.create_group("bins")

            write_bins(res, np_chroms_names,
                       np_chroms_length, grp_bins, h5_opts)

            cell_groups = res_grp.create_group("cells")

            for i in range(cell_cnt):
                print("cell "+str(i)+":")
                cur_cell_grp = cell_groups.create_group(
                    "cell_id"+str(i))
                cell_grp_chroms = cur_cell_grp.create_group("chroms")
                cell_chroms_names, cell_chroms_length = write_chroms(
                    cell_grp_chroms, self.hic, h5_opts)

                cell_grp_bins = cur_cell_grp.create_group("bins")
                write_bins(res, np_chroms_names,
                           np_chroms_length, cell_grp_bins, h5_opts)

                n_chroms = len(cur_cell_grp["chroms"].get("length"))
                n_bins = len(cur_cell_grp["bins"].get("chrom"))

                cell_grp_pixels = cur_cell_grp.create_group("pixels")

                setup_pixels(cell_grp_pixels, n_bins, h5_opts)
                chrom_offset = get_bin_index(
                    cur_cell_grp["bins"], n_chroms, n_bins)

                write_pixels(cell_grp_pixels, self.hic, np_chroms_names, cell_chroms_length,
                             chunksize, res, list(cur_cell_grp["pixels"]), chrom_offset)

                sort_pixel_by_bin(
                    cell_grp_pixels, list(cur_cell_grp["pixels"]))

                n_pixels = len(cur_cell_grp["pixels"].get("bin1_id"))
                bin_offset = get_pixel_index(
                    cur_cell_grp["pixels"], n_bins, n_pixels)
                grp_index = cur_cell_grp.create_group("indexes")

                write_index(grp_index, chrom_offset, bin_offset, h5_opts)


def partition(start, stop, step):
    return ((i, min(i + step, stop)) for i in range(start, stop, step))


class MatrixParser:
    def __init__(self, hic, chrom_names, chrom_lens, chunksize, res, chrom_offset):
        self.hic = hic
        self.chrom_names = chrom_names
        self.chunksize = chunksize
        self.res = res
        self.chrom_offset = chrom_offset
        self.matrix_size = [[j*res*chunksize for j in range(int((chrom_lens[i])/(res*chunksize))+2)]
                            for i in range(chrom_lens.size)]

    def __iter__(self):
        res = self.res
        for index1 in range(len(self.chrom_names)):
            for index2 in range(len(self.chrom_names)):
                print(self.chrom_names[index1], self.chrom_names[index2])
                chrom1 = self.matrix_size[index1]
                chrom2 = self.matrix_size[index2]
                cur_matrix_object = self.hic.getMatrixZoomData(
                    self.chrom_names[index1], self.chrom_names[index2], "observed", "NONE", "BP", self.res)  # KR, SCALE

                for startx in chrom1[:-1]:
                    endx = startx+self.res*self.chunksize
                    for starty in chrom2[:-1]:

                        endy = starty+self.res*self.chunksize

                        cur_numpy_matrix = cur_matrix_object.getRecordsAsMatrix(
                            startx, endx, starty, endy)

                        i, j = np.nonzero(cur_numpy_matrix)

                        yield {
                            "bin1_id": self.chrom_offset[index1]+int(startx/self.res)+i,
                            "bin2_id": self.chrom_offset[index2]+int(starty/self.res)+j,
                            "count": cur_numpy_matrix[i, j],
                        }


if __name__ == "__main__":
    
    root_dir = Path(__file__).resolve().parent.parent.parent.parent
    test = SCHiCGenerator(str(root_dir)+"/hic_data/ENCFF718AWL.hic")
    test.create_all_h5(str(root_dir)+"/hic_data/scHiC5.h5",
                       [100000, 500000, 50000])  # 250000, 100000, 50000, 25000, 10000, 5000
    '''
    root_dir = Path(__file__).resolve().parent.parent.parent.parent
    file_path = str(root_dir)+"/hic_data/scHiC5.h5"

    with h5py.File(file_path, 'r') as hdf:
        all_res = hdf.get("resolutions")
        print(list(all_res))
        all_cell = hdf.get("resolutions/100000/cells")
        print(list(all_cell))
        cur_cell_grp = hdf.get("resolutions/100000/cells/cell_id0")
        n_chroms = cur_cell_grp["chroms"].get("length")
        print(list(n_chroms))

    '''
'''
p1 = []
p2 = []


with h5py.File("data/scHiC.h5", 'r') as hdf:
    G1 = hdf.get("resolutions/100000/cells/cell_id0")
    pixels = G1.get("indexes")

    dset = np.array(pixels.get("bin1_offset"))
    p1 = dset

with h5py.File("data/scHiC2.h5", 'r') as hdf:
    G1 = hdf.get("resolutions/100000/cells/cell_id0")
    pixels = G1.get("pixels")

    dset = np.array(pixels.get("bin1_id"))
    p2 = dset

#print((p1 == p2).all())
# print(p1[:100])
print(p2[:1000])

mzd = hic.getMatrixZoomData('1', '1', "observed", "KR", "BP", 100000)
numpy_matrix = mzd.getRecordsAsMatrix(0, 200000000, 0, 200000000).astype(int)

print(numpy_matrix[0])
i, j = np.nonzero(numpy_matrix)
print(i)
'''
#print(numpy_matrix[i, j].reshape(1, -1)[:100])
