#
# This script performs the GLM analysis on the pre-processed Neuroimaging data.
#
# Useful Links:
#   https://bic-berkeley.github.io/psych-214-fall-2016/classes_and_labs.html
#   https://bic-berkeley.github.io/psych-214-fall-2016/voxel_correlation_solution.html
#   http://practical-neuroimaging.github.io/on_convolution.html
#   http://www.jarrodmillman.com/rcsds/lectures/convolution_background.html
#   https://bic-berkeley.github.io/psych-214-fall-2016/day_07.html
#   http://imaging.mrc-cbu.cam.ac.uk/imaging/Introduction_to_fMRI_2010?action=AttachFile&do=get&target=Intro_fMRI_2010_02_GLM.pdf
#   https://users.fmrib.ox.ac.uk/~stuart/thesis/chapter_6/section6_3.html

import pickle
import nibabel as nib
import numpy as np
import numpy.linalg as npl
from scipy.stats import gamma
from scipy.stats import t as t_dist
import matplotlib.pyplot as plt

class PyGLM(object):
  def __init__(self, bold_data_path, stimulus_file_path, TR, num_vols):
    """
    Args:
      bold_data_path (str): Path to the NIFTI1 data in nii.gz format.
      stimulus_file_path (str): Path to the stimulus file covariates.
      TR (float): Time taken for to record each brain scan.
      num_vols (int): Number of brain volumes acquired during experiment.
    """
    self._ni_data = nib.load(bold_data_path)
    self._bold_data = self._ni_data.get_fdata()
    self._stimulus_file_path = stimulus_file_path
    self._num_vols = num_vols
    self._resolution = TR
    self._contrast = np.atleast_2d([0, 1]).T
    # Total time of motor experiment to get "num_vols" scans of brain where each
    # brain scan was obtained in TR secs.
    self._exp_time_duration = TR * num_vols


  def _get_hrf_response(self, time_duration):
    """
    Links:
      https://practical-neuroimaging.github.io/on_convolution.html
      http://www.jarrodmillman.com/rcsds/lectures/convolution_background.html

    Args:
      time_duration (int): Time duration till which an HRF signal exists.
      step (int): Resolution of the signal in time_duration.

    Returns:
      numpy.ndarray
    """
    time_points = np.arange(0, time_duration, self._resolution)
    hrf_b1 = gamma.pdf(time_points, 6)
    hrf_b2 = gamma.pdf(time_points, 10)
    hrf_signal = hrf_b1 - 0.35*hrf_b2
    return 0.6 * hrf_signal

  def _get_stimulus_signal(self):
    """
    Return a stimulus signal with following specification mentioned in file.
    "Onset" "Duration" "Weight".

    Args:
      stimulus_file_path (str): Covariate file path.

    Returns:
      numpy.ndarray
    """
    f = open(self._stimulus_file_path)
    lines = f.readlines()
    f.close()
    onset = []
    duration = []
    weight = []
    for line in lines:
      line = line.split()
      onset.append(float(line[0]))
      duration.append(float(line[1]))
      weight.append(float(line[2]))

    # Create a time axis with 0.1 sec resolution.
    time_axis = np.arange(0, self._exp_time_duration, self._resolution)
    time_axis[:] = 0.0
    for i in range(len(time_axis)):
      for onst, drtn in zip(onset, duration):
        time_axis[int(onst / self._resolution) :
                  int((onst + drtn) / self._resolution)] = 1

    return time_axis

  def _get_hrf_signal(self):
    """
    Returns the convolved HRF response with Stimulus signal.
    """
    hrf_response = self._get_hrf_response(20)
    stimulus_signal = self._get_stimulus_signal()
    # We need convolved signal only up till self._exp_time_duration. The total
    # length of the convolved signal is len(stimulus_signal) + len(hrf_response)
    # - 1. We need convolved signal up till length len(stimulus_signal).

    # Also note that np.convolve(x,y) = np.convolve(y,x).
    return np.convolve(stimulus_signal, hrf_response)[:-(len(hrf_response) - 1)]

  def _get_design_matrix(self):
    """
    Returns a design matrix X in Y = AX + E, where Y is the bold data, A is the
    coefficient matrix to be estimated and E is the error matrix.

    Returns:
      numpy.ndarray
    """
    # The dimension of the design matrix would be number of brain scans X
    # number of regressors. As of now only two regressors are kept, with one
    # being the intercept and other being the task regressor.
    design_mat = np.zeros((self._num_vols, 2))

    # Set the first column as the intercept.
    design_mat[:,0] = 1
    # Set the second column as the task regressor.
    design_mat[:,1] = self._get_hrf_signal()

    return design_mat

  def _get_error_matrix(self, num_noise_sources=1):
    """
    Generates a random gaussian noise matrix of dimension = number of brain
    scans X number of noise sources. We assume errors to have independent and
    identical normal distributions around 0.

    Args:
      num_noise_sources (int): By default only one noise source is assumed.

    Returns:
      numpy.ndarray
    """
    error_mat = np.zeros((self._num_vols, num_noise_sources))
    for i in range(num_noise_sources):
      error_mat[:,i] = np.random.normal(size=self._num_vols)

    return error_mat

  def _do_glm_for_one_voxel(self, v_x, v_y, v_z):
    """
    Does a GLM analysis for one voxel with coordinate (v_x, v_y, v_z).

    Args:
      v_x (int): x coordinate of the voxel.
      v_y (int): y coordinate of the voxel.
      v_z (int): z coordinate of the voxel.

    Returns:
      TODO
    """
    voxel_time_course = self._bold_data[v_x, v_y, v_z, :] # Y data.
    design_mat = self._get_design_matrix() # X data.
    beta_cap = np.dot(npl.pinv(design_mat), voxel_time_course)
    #beta_cap = np.dot(np.dot(npl.inv(np.dot(design_mat.T, design_mat)),
    #                  design_mat.T), voxel_time_course)
    return beta_cap

  def _get_t_statistic_for_one_voxel(self, v_x, v_y, v_z):
    """
    Returns a T-Statistic score.

    Args:
      v_x (int): x coordinate of the voxel.
      v_y (int): y coordinate of the voxel.
      v_z (int): z coordinate of the voxel.

    Returns:
      TODO
    """
    design_mat = self._get_design_matrix() # X data.
    voxel_time_course = self._bold_data[v_x, v_y, v_z, :] # Y data.
    beta_cap = self._do_glm_for_one_voxel(v_x, v_y, v_z)
    Y_cap = np.dot(design_mat, beta_cap) # Estimated voxel_time_course

    # Error in estimation of voxel's time series.
    error = voxel_time_course - Y_cap
    RSS = (error**2).sum() # Residual Sum of Squares.
    # Degrees of Freedom = Num. of observations - Num. of independent regressors
    # Here "Num. of observations" is simply the number of observations made for
    # a voxel, thus the number of rows in the design matrix (which is equal to
    # the number of rows in voxel_time_course matrix). Next, "Num. of
    # independent regressors" is equal to the rank of design matrix.
    dof = design_mat.shape[0] - npl.matrix_rank(design_mat)
    MRSS = float(RSS) / dof # Mean RSS.
    # Calcuate the Standard Error.
    SE = np.sqrt(MRSS * np.dot(self._contrast.T, np.dot(npl.inv(np.dot(
                 design_mat.T, design_mat)), self._contrast)))
    t_stat = np.dot(self._contrast.T, beta_cap) / SE

    # Get p value for t value using CDF of t distribution.
    lower_tail_p = t_dist.cdf(t_stat, dof)
    upper_tail_p = 1 - lower_tail_p
    return beta_cap, t_stat, dof, upper_tail_p

if __name__=="__main__":
  stimulus_file_path = (
      "/Personal/projects/COL786/data/data_5/subdata/covariates/left_f")
  bold_data_path = (
      "/Personal/projects/COL786/data/data_5/subdata/func/"
      "sub-MSC01_ses-func01_task-motor_run-01_bold.nii.gz")
  processed_bold_data_path = ("/Personal/projects/COL786/Assignment6/"
                              "python_glm_data.feat/filtered_func_data.nii.gz")

  ob = PyGLM(processed_bold_data_path, stimulus_file_path, 2.2, 104)

  v_x, v_y, v_z = 31, 33, 29
  design_mat = ob._get_design_matrix()
  beta_cap, t_stat, dof, upper_tail_p = ob._get_t_statistic_for_one_voxel(
      v_x, v_y, v_z)
  plt.plot(np.dot(design_mat, beta_cap))
  plt.plot(ob._bold_data[v_x, v_y, v_z, :])
  plt.show()
#  t_stat_list = []
#  p_list = []
#  beta_cap_list = []
#  for v_x in range(ob._bold_data.shape[0]):
#    for v_y in range(ob._bold_data.shape[1]):
#      for v_z in range(ob._bold_data.shape[2]):
#        voxel_time_course = ob._bold_data[v_x, v_y, v_z, :]
#        beta_cap, t_stat, dof, upper_tail_p = ob._get_t_statistic_for_one_voxel(
#            v_x, v_y, v_z)
#        #print ("Beta: ", beta_cap, "T-Stat: ", t_stat, "DOF: ", dof, "P: ",
#        #       upper_tail_p)
#        #plt.plot(np.dot(design_mat, beta_cap), label="Assumed HRF signal")
#        #plt.plot(voxel_time_course, label="Voxel's Time Series")
#        #plt.show()
#        t_stat_list.append((t_stat, (v_x, v_y, v_z)))
#        beta_cap_list.append(beta_cap)
#        p_list.append(upper_tail_p)
#    print "x_dim: %s Done!" % v_x
#
#  pickle.dump(t_stat_list, open("T-Stat_list.p", "wb"))
#  pickle.dump(beta_cap_list, open("beta_cap_list.p", "wb"))
#  pickle.dump(p_list, open("p_list.p", "wb"))
