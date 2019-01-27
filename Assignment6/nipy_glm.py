#
# This script used nipy to perform the GLM analysis of the pre-processed
# Neuroimaging data.
#
# Links: https://nipype.readthedocs.io/en/latest/users/examples/fmri_nipy_glm.html#preprocessing-pipeline-nodes
#


import nipy as nip
import numpy as np
import pickle

from nipy.modalities.fmri.glm import GeneralLinearModel
from nipy.modalities.fmri.design import block_amplitudes

class NipyGLM(object):
  def __init__(self, prcsd_bold_data_path, stmls_file_path, TR):
    """
    prcsd_bold_data_path (str): Path to the NIFTI1 file in nii.gz format.
    stmls_file_path (str): Path to the stimulus file covariate.
    TR (float): Time take to record each brain scan.
    """
    self._ni_data = nip.load_image(prcsd_bold_data_path).get_data()
    self._sf_path = stmls_file_path
    self._num_vols = self._ni_data.shape[3]
    self._TR = TR
    # Get the block specification where the 2D data columns represent "start",
    # "end", "amplitude".
    self._block_spec = np.loadtxt(self._sf_path) # start, duration, amplitude.
    # Get "end" from "start" + "duration".
    self._block_spec[0][1] += self._block_spec[0][0]
    self._block_spec[1][1] += self._block_spec[1][0]
    self._task_regressor = block_amplitudes(
        "left_f", self._block_spec, np.arange(self._num_vols) * self._TR)[0]
    self._intercept = np.ones(self._num_vols)
    # Prepare the design matrix: X with dimension: # brain vols x # regressors.
    self._design_matrix = np.vstack((self._intercept, self._task_regressor)).T
    # 1st column in design matrix is an intercept, 2nd column is task regressor.
    self._contrast = np.hstack((0, 1))

  def _do_nipy_glm_for_one_voxel(self, v_x, v_y, v_z):
    """
    Does a NIPY GLM analysis for one voxel with coordinate (v_x, v_y, v_z).

    Args:
      v_x (int): x coordinate of the voxel
      v_y (int): y coordinate of the voxel
      v_z (int): z coordinate of the voxel
    """
    voxel_time_course = np.atleast_2d(self._ni_data[v_x, v_y, v_z, :]).T # Y data.
    model = GeneralLinearModel(self._design_matrix)
    model.fit(voxel_time_course)
    z_val = model.contrast(self._contrast).z_score()
    return z_val

if __name__ == "__main__":
  stimulus_file_path = (
      "/Personal/projects/COL786/data/data_5/subdata/covariates/left_f")
  processed_bold_data_path = ("/Personal/projects/COL786/Assignment6/"
                              "python_glm_data.feat/filtered_func_data.nii.gz")
  ob = NipyGLM(processed_bold_data_path, stimulus_file_path, 2.2)
  #v_coord = (31, 33, 29)
  #ob._do_nipy_glm_for_one_voxel(*v_coord)

  z_stat_list = []
  shape = ob._ni_data.shape # Get the range of x, y, z.
  for v_x in range(shape[0]):
    for v_y in range(shape[1]):
      for v_z in range(shape[2]):
        z_val = ob._do_nipy_glm_for_one_voxel(v_x, v_y, v_z)
        z_stat_list.append((z_val, (v_x, v_y, v_z)))
  pickle.dump(z_stat_list, open("nipy_glm_stat_list.p", "wb"))
