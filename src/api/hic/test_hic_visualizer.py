import time
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path
import h5py

from .hic_visualizer import parse_region, hic_fetch_maps
matplotlib.use('Agg')

REDMAP = LinearSegmentedColormap.from_list(
    "bright_red", [(1, 1, 1), (1, 0, 0)])


def test_visualize_map(fpath, range1, range2):
    # c = Hic_wrapper(cool_uri)

    # future 1: check if use chroms size for end

    row_chrom, row_lo, row_hi = parse_region(range1)
    col_chrom, col_lo, col_hi = parse_region(range2)

    # future 2: check matrix size
    # file_path, root_path, cell_id
    cnt = 1
    cell_ids = [str(i) for i in range(cnt)]

    start_time = time.time()
    arr = hic_fetch_maps(fpath, "100000", cell_ids, range1, range2, True)
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

def test_embeddings(file_path):
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

if __name__ == "__main__":


    # hic_get_chrom_len('scHiC5.h5', '50000', '0')

    root_path = Path(__file__).resolve().parent.parent.parent.parent
    file_path = str(root_path)+"/hic_data/mouse2_slice99.h5"

    test_visualize_map(file_path, "chrom1:0-19515427", "chrom1:0-19515427")