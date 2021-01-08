import numpy as np
import os
import sys
import starfish
from starfish import Experiment
from starfish.types import Features, Axes
from starfish import IntensityTable
from starfish.image import Filter
from starfish.spots import DetectPixels
from datetime import datetime


def DARTFISH_pipeline(fov, codebook, magnitude_threshold, normalize,
                      area_threshold=(5, 100), distance_threshold=3):
    imgs = fov.get_image(starfish.FieldOfView.PRIMARY_IMAGES)

    gauss_filt = Filter.GaussianLowPass(0.7, True)
    gauss_imgs = gauss_filt.run(imgs)

    sc_filt = Filter.Clip(p_max=100, expand_dynamic_range=True)
    norm_imgs = sc_filt.run(gauss_imgs)

    z_filt = Filter.ZeroByChannelMagnitude(thresh=.05, normalize=normalize)
    filtered_imgs = z_filt.run(norm_imgs)

    def compute_magnitudes(stack, norm_order=2):
        pixel_intensities = IntensityTable.from_image_stack(stack)
        feature_traces = pixel_intensities.stack(traces=(Axes.CH.value, Axes.ROUND.value))
        norm = np.linalg.norm(feature_traces.values, ord=norm_order, axis=1)
        return norm

    mags = compute_magnitudes(filtered_imgs)

    psd = DetectPixels.PixelSpotDecoder(
        codebook=codebook,
        metric='euclidean',
        distance_threshold=distance_threshold,
        magnitude_threshold=magnitude_threshold,
        min_area=area_threshold[0],
        max_area=area_threshold[1]
    )

    spot_intensities, results = psd.run(filtered_imgs)
    spot_intensities = IntensityTable(spot_intensities.where(spot_intensities[Features.PASSES_THRESHOLDS], drop=True))
    # reshape the spot intensity table into a RxC barcode vector
    pixel_traces = spot_intensities.stack(traces=(Axes.ROUND.value, Axes.CH.value))

    # extract dataframe from spot intensity table for indexing purposes
    pixel_traces_df = pixel_traces.to_features_dataframe()
    pixel_traces_df['area'] = np.pi * pixel_traces_df.radius ** 2
    return pixel_traces_df, mags


def process_experiment(experiment: starfish.Experiment, output_dir, magnitude_threshold, normalize,
                       area_threshold=(5, 100), distance_threshold=3):
    decoded_intensities = {}
    regions = {}
    count = 0
    for i, (name_, fov) in enumerate(experiment.items()):
        print(datetime.now().strftime(
            '%Y-%d-%m_%H:%M:%S: Started Processing FOV {:03d} with Barcode Magnitude threshold {}'.format(count,
                                                                                                          magnitude_threshold)))
        pixel_traces_df, mags = DARTFISH_pipeline(fov, experiment.codebook, magnitude_threshold, normalize,
                                                  area_threshold, distance_threshold)
        pixel_traces_df.to_csv(
            os.path.join(output_dir, 'starfish_table_bcmag_{}_FOV{:03d}'.format(magnitude_threshold, count) + '.csv'))
        print(datetime.now().strftime(
            '%Y-%d-%m_%H:%M:%S: Finished Processing FOV {:03d} with Barcode Magnitude threshold {}'.format(count,
                                                                                                           magnitude_threshold)))
        count += 1


def starfish_decode(output_dir, decode_dir, magnitude_thresholds=[2.0, 0.9],
                    area_threshold=(5, 100), distance_threshold=3, normalize=True):
    currentTime = datetime.now()
    reportFile = os.path.join(output_dir, currentTime.strftime("%Y-%d-%m_%H:%M_starfish.log"))
    sys.stdout = open(reportFile, 'x')  # redirecting the stdout to the log file
    exp = Experiment.from_json(os.path.join(decode_dir, "experiment.json"))

    for magnitude_threshold in magnitude_thresholds:
        output_path = os.path.join(output_dir, "bcmag{}".format(magnitude_threshold))
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        print(datetime.now().strftime(
            '%Y-%d-%m_%H:%M:%S: Started Processing Experiment with Barcode Magnitude threshold ' + str(
                magnitude_threshold)))
        sys.stdout.flush()
        process_experiment(exp, output_path, magnitude_threshold, normalize=normalize,
                           area_threshold=area_threshold, distance_threshold=distance_threshold)
        print(datetime.now().strftime(
            '%Y-%d-%m_%H:%M:%S: Finished Processing Experiment with Barcode Magnitude threshold ' + str(
                magnitude_threshold)))
        sys.stdout.flush()
    # restoring the stdout pipe to normal
    sys.stdout = sys.__stdout__
