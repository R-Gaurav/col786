#Embedded file name: /Personal/projects/COL786/Assignment4/spatialSmoothing.py
import nibabel as nib
import numpy as np
import sys

def _get_sigma_from_fwhm(fwhm):
  """
  Calculates the corresponding sigma value from `fwhm` value.

  Args:
    fwhm (float): A FWHM value.

  Returns:
    float: The corresponding sigma value.
  """
  return fwhm * 1.0 / np.sqrt(8 * np.log(2))


def _get_finite_gaussian_kernel(x_mean, y_mean, sigma, pix_dim, k_len):
  """
  Args:
    x_mean (float): Mean value of x at which kernel is centered.
    y_mean (float): Mean value of y at which kernel is centered.
    sigma (float): standard deviation of the gaussian kernel.
    pix_dim (float): dimension of the pixel.
    k_len (ing): kernel length (x and y dimensions of the kernel).

  Returns:
    A 2D Gaussian kernel, truncated to 5 * standard deviation. After 3 * standard
    deviation, the values in Gaussian kernel all approach to zero.
  """
  x_vec = np.arange(-(k_len / 2), k_len / 2 + 1, 1) * pix_dim
  y_vec = np.arange(-(k_len / 2), k_len / 2 + 1, 1) * pix_dim
  xx, yy = np.meshgrid(x_vec, y_vec, indexing='ij')
  gauss_2D = np.exp(-(((xx - x_mean) ** 2 + (yy - y_mean) ** 2) / (2 * sigma ** 2)))
  return gauss_2D / np.sum(gauss_2D)


def _do_ifft_fft_convolution(gauss_2d, mat):
  """
  Args:
    gauss_2d (numpy.ndarray): A gaussian kernel.
    mat (numpy.ndarray): A numpy matrix.

  Returns:
    numpy.ndarray: A smoothed matrix.
  """
  print 'Gauss Shape: ', gauss_2d.shape
  print 'mat shape: ', mat.shape
  fr = np.fft.fft2(gauss_2d)
  fr2 = np.fft.fft2(np.flipud(np.fliplr(mat)))
  m, n = fr.shape
  cc = np.real(np.fft.ifft2(fr * fr2))
  cc = np.roll(cc, -m / 2 + 1, axis=0)
  cc = np.roll(cc, -n / 2 + 1, axis=1)
  return cc


def do_gaussian_smoothing(img_file, sigma):
  """
  Does Gaussian smoothing with value of sigma obtained from fwhm.

  Args:
    img_file (str): Path/to/the/NIFTI/*.nii.gz.
    fwhm (float): FWHM parameter.

  Returns:
    numpy.ndarray(:,:,:,:): A 4D array of smoothed value.
  """
  img = nib.load(img_file)
  img_data = img.get_fdata()
  shape = img_data.shape
  ret_mat = np.ndarray(shape)
  gauss_2d = _get_finite_gaussian_kernel(0, 0, sigma, 3.0, 71)
  for vol in range(shape[3]):
      for z in range(shape[2]):
          print 'Value of z: %s, vol: %s' % (z, vol)
          ret_mat[0:71, 0:71, z, vol] = np.round(_do_ifft_fft_convolution(gauss_2d, img_data[0:71, 0:71, z, vol]))
          ret_mat[1:72, 1:72, z, vol] = np.round(_do_ifft_fft_convolution(gauss_2d, img_data[1:72, 1:72, z, vol]))
          ret_mat[1:72, 0:71, z, vol] = np.round(_do_ifft_fft_convolution(gauss_2d, img_data[1:72, 0:71, z, vol]))
          ret_mat[0:71, 1:72, z, vol] = np.round(_do_ifft_fft_convolution(gauss_2d, img_data[0:71, 1:72, z, vol]))

  return ret_mat


if __name__ == '__main__':
  img_file = str(sys.argv[1])
  fwhm = float(sys.argv[2])
  output = str(sys.argv[3])
  sigma = _get_sigma_from_fwhm(fwhm)
  result = do_gaussian_smoothing(img_file, sigma)
  img = nib.Nifti1Image(result, np.eye(4))
  nib.save(img, '/Personal/projects/COL786/data/data_4/sampleResult/' + output + '.nii.gz')
