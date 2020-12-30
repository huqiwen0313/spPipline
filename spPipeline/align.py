import re, shutil, warnings, sys
from time import time
from datetime import datetime
from os import chdir, listdir, getcwd, path, makedirs, remove, walk
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from spPipeline.code_lib import tifffile as tiff # Qiwen: packed as python package, same below
from spPipeline.code_lib import TwoDimensionalAligner_2 as myAligner # Kian: added 201011


def listdirectories(directory='.', pattern='*'):
    """ returns a list of all directories in path"""

    # Qiwen's comment: code is simplified here
    directories = [i for i in next(walk(directory))[1]]
    directories.sort()
    return directories


def mip_gauss_tiled(rnd, fov, dir_root, dir_output='./MIP_gauss',
                        sigma=0.7, channel_int='ch00'):
    """Modified from Matt Cai's MIP.py Maximum intensity projection along z-axis """
    # get current directory and change to working directory - Qiwen's comment don't need this
    # Qiwen: I modified all the dirs so that we can generalize it
    anchor_dir = dir_root + "/" + rnd

    # get all files for position for channel
    image_names = [f for f in listdir(anchor_dir) if
                    re.match(r'.*_s' + '{:02d}'.format(fov) + r'.*_' + channel_int + r'\.tif', f)]

    # put images of correct z_range in list of array
    nImages = len(image_names)
    image_list = [None] * nImages
    for i in range(len(image_names)):
        # Qiwen: modify path here to find target image
        image_list[i] = (ndimage.gaussian_filter(plt.imread(path.join(anchor_dir, image_names[i])),
                                                 sigma=sigma))
    image_stack = np.dstack(image_list)
    max_array = np.amax(image_stack, axis=2)

    # PIL unable to save uint16 tif file
    # Need to use alternative (like libtiff)
    # Make directories if necessary
    if not dir_output.endswith('/'):
        dir_output = "{0}/".format(dir_output)

    if not path.isdir(dir_output + 'FOV{:03d}'.format(fov)):
        makedirs(dir_output + 'FOV{:03d}'.format(fov))

    tiff.imsave(path.join(dir_output, 'FOV{:03d}'.format(fov) + '/MIP_' + rnd + '_FOV{:03d}'.format(fov) +
                    '_' + channel_int + '.tif'), max_array)


class ImageAlign:
    """Class for image alignment and registration """

    def __init__(self, raw_dir, output_dir, rnd_list, n_fovs, sigma,
                 channel_DIC_reference, channel_DIC, cycle_other, channel_DIC_other, **kwargs):
        # path to raw files (e.g. ./0_Raw)
        self.raw_dir = raw_dir
        # path of output files
        self.output_dir = output_dir
        # rounds
        self.rnd_list = rnd_list
        # Number of FOVs (field of views)
        self.n_fovs = n_fovs
        # signma parameter for gaussian filter
        self.sigma = sigma
        self.channel_DIC_reference = channel_DIC_reference
        self.channel_DIC = channel_DIC
        self.cycle_other = cycle_other
        self.channel_DIC_reference = channel_DIC_reference
        self.channel_DIC_other = channel_DIC_other
        self.dir_output_Projected = path.join(output_dir, "1_Projected")
        self.dir_output_aligned = path.join(output_dir, "2_Registered")
        self.cycle_reference = rnd_list[round(len(self.rnd_list) / 2)]

    def get_maximum_intensity(self):
        # MIP
        for rnd in self.rnd_list:
            if "DRAQ5" in rnd or "anchor" in rnd:
                channel_list = [0, 1]
            else:
                channel_list = [0, 1, 2, 3]

            for channel in channel_list:
                print('Generating MIPs for ' + rnd + ' channel {0} ...'.format(channel))
                for fov in range(self.n_fovs):
                    mip_gauss_tiled(rnd, fov, self.raw_dir, self.dir_output_Projected,
                                    sigma=self.sigma, channel_int="ch0{0}".format(channel))
                print('Done\n')

    def dimension_align_2d(self):
        position_list = listdirectories(path.join(self.dir_output_Projected))
        if not path.isdir(self.dir_output_aligned):
            makedirs(self.dir_output_aligned)

        currentTime = datetime.now()
        reportFile = path.join(self.dir_output_aligned,
                               currentTime.strftime("%Y-%d-%m_%H:%M_SITKAlignment.log"))
        # redirecting the stdout to the log file
        sys.stdout = open(reportFile, 'w')

        for position in position_list:
            for rnd in self.rnd_list:
                print(datetime.now().strftime("%Y-%d-%m_%H:%M:%S: " + str(position) +
                                              ', cycle ' + rnd + ' started to align'))
                aligner = myAligner.TwoDimensionalAligner(
                    destinationImagesFolder=path.join(self.dir_output_Projected, position),
                    originImagesFolder=path.join(self.dir_output_Projected, position),
                    originMatchingChannel=self.channel_DIC if rnd not in self.cycle_other else self.channel_DIC_other[rnd],
                    destinationMatchingChannel=self.channel_DIC_reference,
                    imagesPosition=position,
                    destinationCycle=self.cycle_reference,
                    originCycle=rnd,
                    resultDirectory=path.join(self.dir_output_aligned, position),
                    MaximumNumberOfIterations=400)

                for file in [file for file in listdir() if file.startswith('IterationInfo.0')]:
                    if path.isfile(path.join(self.dir_output_aligned, position, "MetaData", file)):
                        remove(path.join(self.dir_output_aligned, position, "MetaData", file))
                        shutil.move(src=file, dst=path.join(self.dir_output_aligned, position, "MetaData"))
