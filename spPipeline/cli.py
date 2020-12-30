import os
import re
import argparse
from starfish.types import Axes
from spPipeline.codebookGenerator import create_jsoncodebook
from spPipeline.align import *
from spPipeline.StitchDriver import Stitch
from spPipeline.toStarfishFormat import format_data
from spPipeline.starfishDecode import *
from spPipeline.combineFOVs import combine_fovs
from spPipeline.segmentation import segmentation

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--raw", default="../0_Raw", help="input dir contains Raw files")
parser.add_argument("-o", "--output", default="./", help="output dir")
parser.add_argument("-s", "--sigma", default=0.7, type=float,
                    help="sigma for gaussian filter")
parser.add_argument("-rd", "--rnd_list",
                    default=["0_anchor", "1_dc0", "2_dc1", "3_dc2", "4_dc3", "5_dc4", "6_dc5", "7_DRAQ5"],
                    type=list,
                    help="list of decoding cycle rounds")
parser.add_argument("-nf", "--nfovs", default=80, type=int, help="number of fov (filed of view)")
parser.add_argument("-cr", "--channel_DIC_reference", default='ch03',
                    help="DIC channel for reference cycle")
parser.add_argument("-cd", "--channel_DIC", default='ch03',
                    help="DIC channel for (non-reference) decoding cycles")
parser.add_argument("-co", "--cycle_other", type=list, default=['0_anchor', '7_DRAQ5'],
                    help="other data-containing folders which need to be aligned but are not names 'CycleXX'")
parser.add_argument("-cdo", "--channel_DIC_other", default={'0_anchor': 'ch01', '7_DRAQ5': 'ch01'},
                    help="DIC channel for other data-containing folders")
parser.add_argument("-gx", "--grid_size_x", type=int, default=8,
                    help="grid size (x-axis) for each tile")
parser.add_argument("-gy", "--grid_size_y", type=int, default=10,
                    help="grid size (y-axis) for each tile")
parser.add_argument("-to", "--tile_overlap", type=int, default=15,
                    help="overlap for each tile")
parser.add_argument("-ij", "--ij_path", default="/home/qiwenhu/software/Fiji.app/ImageJ-linux64",
                    help="The path to imagej’s executables")
parser.add_argument("-sr", "--stitchRef", default="dc3",
                    help="The round to be used as the reference for stitching")

args = parser.parse_args()


def main():
    # Synthesizing codebook
    cycle_reference = args.rnd_list[round(len(args.rnd_list) / 2)]
    barcode_file = os.path.join(args.output, '_codebook/TB12k_Mar2018_V7_noAnchor.txt')
    codebook_file = os.path.join(args.output, '_codebook/TB12k_Mar2018_V7_noAnchor.json')
    create_jsoncodebook(infilepath=barcode_file, outfilepath=codebook_file,
                        totalCycles=6, offCycles=2, firstCycleAnchor=False,
                        addEmptyBarcodes=True, uniColorAllowed=False)

    # image align and maximum projection
    image_align = ImageAlign(raw_dir=args.raw, output_dir=args.output,
                             rnd_list=args.rnd_list, n_fovs=args.nfovs, sigma=args.sigma,
                             channel_DIC_reference=args.channel_DIC_reference,
                             channel_DIC=args.channel_DIC, cycle_other=args.cycle_other,
                             channel_DIC_other=args.channel_DIC_other)
    image_align.get_maximum_intensity()
    image_align.dimension_align_2d()

    # stitching
    input_dir = os.path.join(args.output, "2_Registered")
    stitch_dir = os.path.join(args.output, "2_Registered/stitched")
    rounds = [re.sub(r'\d_', r"", i) for i in args.rnd_list]
    stitchChRef = args.channel_DIC_reference
    image_stitching = Stitch(input_dir=input_dir, stitch_dir=stitch_dir, rounds=rounds,
                             stitchRef=args.stitchRef, stitchChRef=stitchChRef,
                             grid_size_x=args.grid_size_x, grid_size_y=args.grid_size_y,
                             tileOverlap=args.tile_overlap, ij_path=args.ij_path)
    image_stitching.stitch_reference()
    image_stitching.stitch_tileconfig()
    image_stitching.generate_cvs()

    # converting to starfish format
    SHAPE = {Axes.Y: 1024, Axes.X: 1024}
    VOXEL = {"Y": 0.144, "X": 0.144, "Z": 0.420}

    RND_LIST = args.rnd_list[1:]
    RND_ALIGNED = RND_LIST[round(len(RND_LIST) / 2) - 1]
    RND_DRAQ5 = RND_LIST[-1]

    decode_dir = os.path.join(args.output, "3_Decoded/data_Starfish")
    codebook_path = os.path.join(args.output, "_codebook", "TB12k_Mar2018_V7_noAnchor.json")
    if not os.path.exists(codebook_path):
        raise FileNotFoundError("Codebook Not Found.")
    if not os.path.exists(decode_dir):
        os.makedirs(decode_dir)
    format_data(input_dir, decode_dir, fov_count=args.nfovs, SHAPE=SHAPE,
                voxel=VOXEL, rnd_list=RND_LIST, rnd_aligned=RND_ALIGNED, rnd_draq5=RND_DRAQ5,
                codebook_path=codebook_path, rounds=6, channels=3, zplanes=1)

    # starfish decoding
    output_dir = os.path.join(args.output, "3_Decoded/output_Starfish")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    starfish_decode(output_dir=output_dir, decode_dir=decode_dir,
                    magnitude_thresholds=[2.0], area_threshold=(5, 100),
                    distance_threshold=3, normalize=True)
    starfish_decode(output_dir=output_dir, decode_dir=decode_dir,
                    magnitude_thresholds=[0.9], area_threshold=(5, 100),
                    distance_threshold=3, normalize=False)

    # Pooling rolonies from all FOVs and filtering
    combine_fovs(decoding_dir=output_dir, voxel=VOXEL, emptyFractionThresh=0.12)

    # cell segmentation
    nuc_path = os.path.join(args.output, "2_Registered/stitched/MIP_7_DRAQ5_ch00.tif")
    saving_path = os.path.join(args.output, '4_CellAssignment')
    bcmag = 'bcmag2.0'
    spot_file = os.path.join(args.output, '3_Decoded/output_Starfish/{}/all_spots_filtered.tsv'.format(bcmag))
    segmentation(nuc_path, saving_path, bcmag, spot_file)
