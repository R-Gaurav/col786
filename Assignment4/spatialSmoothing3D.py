#
# This file implements a 3D smoothing of a volume of brain.
#

import nibabel as nib
import numpy as np
import sys

class SpatialSmoothing(object):
  def __init__(self, fwhm, data_path, out_file_name, xpix_dim, ypix_dim,
               zpix_dim):
    """
    Args:
      fwhm (float): Full Width Half Maximum of the Gaussian Kernel.
      data_path (string): Path/to/data.nii.gz
      out_file_name (string): Output file's name.
    """
    self._fwhm = fwhm
    self._data = nib.load(data_path).get_fdata()
    self._output = out_file_name
    self._sigma = self._get_sigma_from_fwhm(self._fwhm)
    self._xlmt, self._ylmt, self._zlmt, self._nsple = self._data.shape
    self._xpix_dim = xpix_dim
    self._ypix_dim = ypix_dim
    self._zpix_dim = zpix_dim

  def _get_sigma_from_fwhm(self, fwhm):
    """
    Calculates the corresponding sigma value from `fwhm` value.

    Args:
      fwhm (float): A FWHM value.

    Returns:
      float: The corresponding sigma value.
    """
    return self._fwhm * 1.0 / np.sqrt(8 * np.log(2))

  def _get_finite_gaussian_kernel(self, x_mean, y_mean, z_mean):
    """
    Args:
      x_mean (float): Mean value of x at which kernel is centered.
      y_mean (float): Mean value of y at which kernel is centered.
      z_mean (float): Mean value of z at which kernel is centered.

    Returns:
      numpy.ndarray (:,:,:): A 3D Gaussian kernel
    """
    # Get odd values for a central point in a 3D cuboid.
    x_len = self._xlmt - 1
    y_len = self._ylmt - 1
    z_len = self._zlmt - 1

    x_vec = np.arange(-(x_len / 2), x_len / 2 + 1, 1) * self._xpix_dim
    y_vec = np.arange(-(y_len / 2), y_len / 2 + 1, 1) * self._ypix_dim
    z_vec = np.arange(-(z_len / 2), z_len / 2 + 1, 1) * self._zpix_dim

    xx, yy, zz = np.meshgrid(x_vec, y_vec, z_vec, indexing='ij')
    gauss_3D = np.exp(
        -(((xx - x_mean) ** 2 + (yy - y_mean) ** 2 + (zz - z_mean) ** 2) /
         (2 * self._sigma ** 2)))
    return gauss_3D / np.sum(gauss_3D)

  def _do_fft_ifft_convolution(self, gauss_3D, cuboid):
    """
    Args:
      gauss_3D (numpy.ndarray): A 3D Gaussian Kernel.
      cuboid (numpy.ndarray): A 3D brain volume.

    Returns:
      numpy.ndarray (:,:,:): A 3D smoothed array.
    """
    fft1 = np.fft.fftn(gauss_3D)
    fft2 = np.fft.fftn(cuboid)
    m, n, o = fft2.shape

    ifft = np.real(np.fft.ifftn(fft1 * fft2))
    ifft = np.roll(ifft, -m / 2 + 1, axis=0)
    ifft = np.roll(ifft, -n / 2 + 1, axis=1)
    ifft = np.roll(ifft, -o / 2 + 1, axis=2)

    return ifft

  def do_gaussian_smoothing(self):
    """
    Does Gaussian smoothing with value of sigma obtained from fwhm.

    Returns:
      numpy.ndarray(:,:,:,:): A 4D array of smoothed value.
    """
    ret_vol = np.ndarray((self._xlmt, self._ylmt, self._zlmt, self._nsple))
    gauss_3D = self._get_finite_gaussian_kernel(0,0,0)
    for sample in range(self._nsple):
      brain_vol_000 = self._data[0:self._xlmt-1, 0:self._ylmt-1, 0:self._zlmt-1,
                                 sample]
      brain_vol_010 = self._data[0:self._xlmt-1, 1:self._ylmt, 0:self._zlmt-1,
                                 sample]
      brain_vol_100 = self._data[1:self._xlmt, 0:self._ylmt-1, 0:self._zlmt-1,
                                 sample]
      brain_vol_110 = self._data[1:self._xlmt, 1:self._ylmt, 0:self._zlmt-1,
                                 sample]


      brain_vol_001 = self._data[0:self._xlmt-1, 0:self._ylmt-1, 1:self._zlmt,
                                 sample]
      brain_vol_011 = self._data[0:self._xlmt-1, 1:self._ylmt, 1:self._zlmt,
                                 sample]
      brain_vol_101 = self._data[1:self._xlmt, 0:self._ylmt-1, 1:self._zlmt,
                                 sample]
      brain_vol_111 = self._data[1:self._xlmt, 1:self._ylmt, 1:self._zlmt,
                                 sample]


      ret_vol[0:self._xlmt-1, 0:self._ylmt-1, 0:self._zlmt-1, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_000))
      ret_vol[0:self._xlmt-1, 1:self._ylmt, 0:self._zlmt-1, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_010))
      ret_vol[1:self._xlmt, 0:self._ylmt-1, 0:self._zlmt-1, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_100))
      ret_vol[1:self._xlmt, 1:self._ylmt, 0:self._zlmt-1, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_110))


      ret_vol[0:self._xlmt-1, 0:self._ylmt-1, 1:self._zlmt, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_001))
      ret_vol[0:self._xlmt-1, 1:self._ylmt, 1:self._zlmt, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_011))
      ret_vol[1:self._xlmt, 0:self._ylmt-1, 1:self._zlmt, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_101))
      ret_vol[1:self._xlmt, 1:self._ylmt, 1:self._zlmt, sample] = np.round(
          self._do_fft_ifft_convolution(gauss_3D, brain_vol_111))

    return ret_vol

if __name__ == "__main__":
  img_file = str(sys.argv[1])
  fwhm = float(sys.argv[2])
  output = str(sys.argv[3])
  ob = SpatialSmoothing(fwhm, img_file, output, 3.0, 3.0, 3.33)
  result = ob.do_gaussian_smoothing()
  img = nib.Nifti1Image(result, np.eye(4))
  nib.save(img, '/Personal/projects/COL786/data/data_4/sampleResult/3D' +
      output + '.nii.gz')
