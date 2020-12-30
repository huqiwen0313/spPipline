import os, re
import shutil, sys
import pandas as pd
import spPipeline.code_lib.IJ_stitch_201020 as IJS
from datetime import datetime


def add2dict2dict(key, value, dic):
    if key in dic:
        dic[key].append(value)
    else:
        dic[key] = [value]


def copy2dir(files2copy, dest_dir):
    for infile in files2copy:
        shutil.copy2(infile, dest_dir)


def changeTileConfig(reffile, nrefile, nrefNames, fov_pat):
    """ Looping through all lines of the reference tile config, and change the channel
        to match the non-reference image files.

        Args:
            reffile: path to the reference tile configuration file
            nrefile: path to the tile configation file that we want to generate
            nrefNames: file names that need to be substituted for the original reference filenames.
            fov_pat: regex pattern that specifies the FOV.
    """
    with open(nrefile, 'w') as writer, open(reffile, 'r') as reader:
        for line in reader:
            refmtch = re.search(".tif", line)  # assuming all images are .tif
            if refmtch is None:
                writer.writelines(line)
            else:
                fov = re.search(fov_pat, line).group(0)  # the FOV in this line

                # finding the non-ref image with the same fov
                for nrefn in nrefNames:
                    if fov in nrefn:
                        # substituting the whole file name section
                        writer.writelines(re.sub(r"^\S+.tif", nrefn, line))


def cleanUpImages(file_dict, file_dir):
    """ Deleted the images we moved for stitching """
    for key in file_dict:
        for file in file_dict[key]:
            os.remove(os.path.join(file_dir, os.path.basename(file)))


def writeReport(spOut):
    dtn = datetime.now()
    dtn = "{0}-{1}-{2}_{3}:{4}:{5}".format(dtn.year, dtn.month, dtn.day,
                                           dtn.hour, dtn.minute, dtn.second)
    print("{0}: ImageJ's stdout:".format(dtn))
    print(spOut.stdout)
    print("{0}: ImageJ's stderr:".format(dtn))
    print(spOut.stderr)


def readStitchInfo(infoFile, rgx):
    """ read ImageJ's stitching output and spit out the top left position of
    each tile image on the stitched image in a dataframe"""
    with open(infoFile, 'r+') as reader:
        infoDict = {}
        for line in reader:
            if line.startswith('# Define the image coordinates'):
                break

        positions, xs, ys = [], [], []
        for line in reader:
            pos_re = re.search(rgx, line)
            positions.append(pos_re.group('fov'))

            coord_re = re.search(r".tif.*\(([-+]?[0-9]*[.][0-9]*)" +
                                 r".*?([-+]*[0-9]*[.][0-9]*)\)", line)

            xs.append(float(coord_re.group(1)))
            ys.append(float(coord_re.group(2)))

        return pd.DataFrame({'fov': positions, 'x': xs, 'y': ys})


class Stitch:
    """ We want to stitch all channels of all cycles of DART-FISH.
        Since all the images that need to be stitched have to in the same directory,
        we have to move images of different FOVs in the same directory and run the image stitching.
        As of now (Oct 14th, 2020), after maximum projecting and registering, images of the same FOV
        are kept in the same directory.
        This code assumes that all images are registered, so one specified cycle and channel is
        used to find the tile configuration and that setting will be applied to all other tiles and channels.
        IMPORTANT: The ImageJ path has to be set with in arguments
    """
    def __init__(self, input_dir, stitch_dir, rounds, stitchRef, stitchChRef,
                 grid_size_x, grid_size_y, tileOverlap, ij_path):
        self.input_dir = input_dir
        self.stitch_dir = stitch_dir
        self.rounds = rounds
        self.stitchRef = stitchRef
        self.grid_size_x = grid_size_x
        self.grid_size_y = grid_size_y
        self.tileOverlap = tileOverlap
        self.stitchChRef = stitchChRef
        self.ij_path = ij_path

        # look for file patterns
        # 0: all, 1: MIP_rnd#, 2:dc/DRAQ, 3: FOV, 4: chfile_regex = re.compile(filePattern)
        self.filePattern = r"(?P<intro>\S+)?_(?P<rndName>\S+)_(?P<fov>FOV\d+)_(?P<ch>ch\d+)\S*.tif$"
        self.file_regex = re.compile(self.filePattern)
        # pattern to extract the fov number
        self.fov_pat = r"(FOV)\d+"
        # string to substitite the fov# with {iii}
        self.fov_sub = r"\1{iii}"
        self.orig_stdout = sys.stdout

        # get fovs information
        fovs = [file for file in os.listdir(self.input_dir) if re.match("FOV\d+", file)]
        self.fovs = sorted(fovs, key=lambda x: int(x[3:]))

    def stitch_reference(self):
        """Stitch the reference channel with random-image fusion"""

        if not os.path.isdir(self.stitch_dir):
            os.mkdir(self.stitch_dir)

        # Redirecting stdout to write in a report file
        dtn = datetime.now()
        reportfile = os.path.join(self.stitch_dir,
                                  "{0}-{1}-{2}_{3}:{4}:{5}-stitch report.txt".format(dtn.year,
                                                                                     dtn.month, dtn.day,
                                                                                     dtn.hour, dtn.minute, dtn.second))
        reporter = open(reportfile, 'w')
        sys.stdout = reporter

        if not self.stitchRef in self.rounds:
            raise ValueError("Stitching reference round is not in rounds list: {}".format(self.rounds))

        # Copy images to stitching folder
        refImgPaths = []
        for fov in self.fovs:
            fov_files = os.listdir(os.path.join(self.input_dir, fov))

            for file in fov_files:
                mtch = self.file_regex.match(file)
                if mtch is not None:
                    if (mtch.group('rndName') == self.stitchRef) and (mtch.group('ch') == self.stitchChRef):
                        refImgPaths.append(os.path.join(self.input_dir, fov, file))
        copy2dir(refImgPaths, self.stitch_dir)

        print(datetime.now().strftime("%Y-%m-%d_%H:%M:%S: Stitching reference {0}, {1}".format(self.stitchRef, self.stitchChRef)))
        refTileConfigFile = "Ref_{0}_{1}_TileConfig.txt".format(self.stitchRef, self.stitchChRef)
        f_pat = re.sub(self.fov_pat, self.fov_sub, os.path.basename(refImgPaths[0]))  # ImageJ sequence pattern

        refStitcher = IJS.IJ_Stitch(input_dir=self.stitch_dir, output_dir=self.stitch_dir, file_names=f_pat,
                                    imagej_path=self.ij_path, Type='Grid: row-by-row', Order='Left & Up',
                                    tile_overlap=self.tileOverlap, grid_size_x=self.grid_size_x,
                                    grid_size_y=self.grid_size_y,
                                    output_textfile_name=refTileConfigFile,
                                    fusion_method='Intensity of random input tile',
                                    compute_overlap=True,
                                    macroName='{0}_{1}.ijm'.format(self.stitchRef, self.stitchChRef),
                                    output_name='Ref_{0}_{1}_random_fusion.tif'.format(self.stitchRef, self.stitchChRef))
        res = refStitcher.run()
        writeReport(res)

    def stitch_tileconfig(self):
        # Stitch everything using the reference TileConfig
        for rnd in self.rounds:
            # Copy images to stitching folder
            # contains the path to images-to-be-stitched in each round
            thisRnd = {}
            for fov in self.fovs:
                fov_files = os.listdir(os.path.join(self.input_dir, fov))

                for file in fov_files:
                    mtch = self.file_regex.match(file)
                    if mtch is not None:
                        if mtch.group('rndName') == rnd:
                            add2dict2dict(mtch.group('ch'),
                                          os.path.join(self.input_dir, fov, file), thisRnd)

            chans = list(thisRnd)
            for ch in chans:
                copy2dir(thisRnd[ch], self.stitch_dir)

            for nch in chans:
                nrefTileConfig = "{0}-to-{1}_{2}_TileConfig.registered.txt".format(self.stitchRef, rnd, nch)
                changeTileConfig(reffile=os.path.join(self.stitch_dir,
                                                      "Ref_{0}_{1}_TileConfig.registered.txt".format(self.stitchRef,
                                                                                                     self.stitchChRef)),
                                 nrefile=os.path.join(self.stitch_dir, nrefTileConfig),
                                 nrefNames=[os.path.basename(f) for f in thisRnd[nch]],
                                 fov_pat=self.fov_pat
                                 )

                f_pat = re.sub(self.fov_pat, self.fov_sub, os.path.basename(thisRnd[nch][0]))  # ImageJ sequence pattern

                print(datetime.now().strftime(
                    "%Y-%m-%d_%H:%M:%S: Stitching round {0}, {1} using the coordinates from {2}".format(rnd, nch,
                                                                                                        self.stitchRef)))
                nonRefStitcher = IJS.IJ_Stitch(input_dir=self.stitch_dir, output_dir=self.stitch_dir, file_names=f_pat,
                                               imagej_path=self.ij_path, Type='Positions from file',
                                               Order='Defined by TileConfiguration',
                                               layout_file=os.path.join(nrefTileConfig),
                                               compute_overlap=False, macroName='{0}_{1}.ijm'.format(rnd, nch),
                                               fusion_method='Max. Intensity')
                res = nonRefStitcher.run()
                writeReport(res)
            cleanUpImages(thisRnd, self.stitch_dir)

    def generate_cvs(self):
        # Writing a CSV file for the coordinates of the registration reference cycle
        allfiles = os.listdir(self.stitch_dir)
        ref_ch = self.stitchChRef  # if rnd in stitchChRefAlt else stitchChRef
        regRef_tileconfig_file = [f for f in allfiles
                                  if f == "Ref_{0}_{1}_TileConfig.registered.txt".format(self.stitchRef, self.stitchChRef)]
        coords = readStitchInfo(os.path.join(self.stitch_dir, regRef_tileconfig_file[0]), self.filePattern[0:-1])
        coords.to_csv(os.path.join(self.stitch_dir, 'registration_reference_coordinates.csv'), index=False)

        sys.stdout = self.orig_stdout  # restoring the stdout pipe to normal








