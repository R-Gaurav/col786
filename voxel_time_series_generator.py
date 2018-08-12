#
# Voxel Time Series Generator.
#
# Accepts <image_name> <x_coordinate> <y_coordinate> <z_coordinate> and creates
# a file named "output.txt" containing the time series of the voxel denoted by
# the passed coordinates in the image.
#
# Note: Image should be in *.nii.gz format or *.hdr or *.img format.
#
# Run as:
# python voxel_time_series_generator.py <image_name> <x> <y> <z>
#

import nibabel as nib
import sys

def generate_time_series(image_name, x_c, y_c, z_c):
  """
  Generates the time series of voxel at (x_c, y_c, z_c) in image "image_name".
  Saves the time series in "output.txt" file.

  Args:
    image_name (str): Image name in string.
    x_c (int): x coordinate.
    y_c (int): y coordinate.
    z_c (int): z coordinate.
  """
  image = nib.load(image_name) # This image is 4D: (72, 72, 38, 120)
  # Last dimension: 120 is the temporal dimension, that is: 120 scannings are
  # done throughout the time of acquiring the image.

  # Assumption with this data: The orientation and reference space for all the
  #                            3D scans of brain are same.
  image_data = image.get_fdata() # This returns an <numpy.ndarray> of 4 dims.

  # Iterate over the 120 3D volumes of brain to get the time series vector of
  # voxel denoted by (x_c, y_c, z_c).
  voxel_time_series = []

  image_shape = image.get_shape() # (72, 72, 38, 120)
  for t in range(image_shape[3]):
    image3D = image_data[:,:,:,t]
    voxel_time_series.append(image3D[x_c, y_c, z_c])

  f = open("Output.txt", "w")
  for element in voxel_time_series:
    f.write(str(element)+" ")
  f.close()

if __name__ == "__main__":
  image_name = sys.argv[1]
  x_c = int(sys.argv[2])
  y_c = int(sys.argv[3])
  z_c = int(sys.argv[4])

  generate_time_series(image_name, x_c, y_c, z_c)
