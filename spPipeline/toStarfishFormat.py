import io
import json
import os
import zipfile
from typing import Mapping, Tuple, Union

import numpy as np
import requests
from skimage.io import imread
from slicedimage import ImageFormat
import pandas as pd

from starfish import Codebook
from starfish.experiment.builder import FetchedTile, TileFetcher
from starfish.experiment.builder import write_experiment_json
from starfish.types import Axes, Coordinates, Features, Number


from shutil import copy2


class DARTFISHTile(FetchedTile):
	def __init__(self, file_path, shape, voxel):
		self.file_path = file_path
		self.VOXEL = voxel
		self.SHAPE = shape

	@property
	def shape(self) -> Tuple[int, ...]:
		return self.SHAPE

	@property
	def coordinates(self) -> Mapping[Union[str, Coordinates], Union[Number, Tuple[Number, Number]]]:
		fov_dir = os.path.dirname(self.file_path)
		registered_dir = os.path.dirname(fov_dir)
		fov = os.path.basename(fov_dir)

		#read coordinates file
		coordinatesTablePath = os.path.join(registered_dir, "stitched",
											 "registration_reference_coordinates.csv")
		
		if os.path.exists(coordinatesTablePath):
			coordinatesTable = pd.read_csv(coordinatesTablePath)
			if coordinatesTable.x.min() < 0:
				coordinatesTable.x = coordinatesTable.x.subtract(coordinatesTable.x.min())
			if coordinatesTable.y.min() < 0:
				coordinatesTable.y = coordinatesTable.y.subtract(coordinatesTable.y.min())
		
			#find coordinates
			locs= coordinatesTable.loc[coordinatesTable.fov == fov].reset_index(drop=True)

			locs = {
				Coordinates.X: (locs.x[0]*self.VOXEL["X"], (locs.x[0] + self.SHAPE[Axes.X])*self.VOXEL["X"]),
				Coordinates.Y: (locs.y[0]*self.VOXEL["Y"], (locs.y[0] + self.SHAPE[Axes.Y])*self.VOXEL["Y"]),
				Coordinates.Z: (0.0, 10.0),
			}
		else:
			print("Coordinate file did not exist at: {}".format(coordinatesTablePath))
			locs = {
				Coordinates.X: (0.0, 0.001),
				Coordinates.Y: (0.0, 0.001),
				Coordinates.Z: (0.0, 0.001),
			}

		return locs

	def tile_data(self) -> np.ndarray:
		return imread(self.file_path)


class DARTFISHPrimaryTileFetcher(TileFetcher):
	def __init__(self, input_dir, rnd_list, shape, voxel):
		self.input_dir = input_dir
		self.RND_LIST = rnd_list
		self.shape = shape
		self.voxel = voxel

	@property
	def ch_dict(self):
		ch_dict = {0: 'ch00', 1: 'ch02', 2: 'ch01'}
		return ch_dict

	@property
	def round_dict(self):
		round_str = self.RND_LIST
		round_dict = dict(enumerate(round_str))
		return round_dict

	def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
		filename = "MIP_{}_FOV{:03d}_{}.tif".format(self.round_dict[r],
												fov, self.ch_dict[ch])
		file_path = os.path.join(self.input_dir, "FOV{:03d}".format(fov), filename)
		return DARTFISHTile(file_path, shape=self.shape, voxel=self.voxel)


class DARTFISHnucleiTileFetcher(TileFetcher):
	def __init__(self, path, rnd_draq5, shape, voxel):
		self.path = path
		self.RND_DRAQ5 = rnd_draq5
		self.shape = shape
		self.voxel = voxel

	def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
		file_path = os.path.join(self.path, "FOV{:03d}".format(fov),
								 "MIP_{}_FOV{:03d}_ch00.tif".format(self.RND_DRAQ5, fov))
		return DARTFISHTile(file_path, shape=self.shape, voxel=self.voxel)


class DARTFISHbrightfieldTileFetcher(TileFetcher):
	def __init__(self, path, rnd_aligned, shape, voxel):
		self.path = path
		self.RND_ALIGNED = rnd_aligned
		self.shape = shape
		self.voxel = voxel

	def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
		file_path = os.path.join(self.path, "FOV{:03d}".format(fov),
								 "MIP_{}_FOV{:03d}_ch03.tif".format(self.RND_ALIGNED, fov))
		return DARTFISHTile(file_path, shape=self.shape, voxel=self.voxel)


def download(input_dir, url):
	print("Downloading data ...")
	r = requests.get(url)
	z = zipfile.ZipFile(io.BytesIO(r.content))
	z.extractall(input_dir)


def write_json(res, output_path):
	json_doc = json.dumps(res, indent=4)
	print(json_doc)
	print("Writing to: {}".format(output_path))
	with open(output_path, "w") as outfile:
		json.dump(res, outfile, indent=4)


def format_data(input_dir, output_dir, fov_count, SHAPE, voxel, rnd_list, rnd_aligned, rnd_draq5,
				codebook_path, rounds=6, channels=3, zplanes=54):
	if not input_dir.endswith("/"):
		input_dir += "/"

	if not output_dir.endswith("/"):
		output_dir += "/"

	def overwrite_codebook(codebook_path, output_dir):
		copy2(codebook_path, os.path.join(output_dir, "codebook.json"))
		
	# the magic numbers here are just for the ISS example data set.
	write_experiment_json(
		output_dir,
		fov_count,
		ImageFormat.TIFF,
		primary_image_dimensions={
			Axes.ROUND: rounds,
			Axes.CH: channels,
			Axes.ZPLANE: zplanes,
		},
		aux_name_to_dimensions={
			'nuclei': {
				Axes.ROUND: 1,
				Axes.CH: 1,
				Axes.ZPLANE: zplanes,
			},
			'dic': {
				Axes.ROUND: 1,
				Axes.CH: 1,
				Axes.ZPLANE: zplanes,
			},
		},
		primary_tile_fetcher=DARTFISHPrimaryTileFetcher(input_dir, rnd_list=rnd_list, shape=SHAPE, voxel=voxel),
		aux_tile_fetcher={
			"nuclei": DARTFISHnucleiTileFetcher(os.path.join(input_dir), rnd_draq5=rnd_draq5, shape=SHAPE, voxel=voxel),
			"dic": DARTFISHbrightfieldTileFetcher(os.path.join(input_dir), rnd_aligned=rnd_aligned, shape=SHAPE, voxel=voxel)
		},
		# postprocess_func=add_codebook,
		default_shape=SHAPE
	)
	overwrite_codebook(codebook_path, output_dir)

