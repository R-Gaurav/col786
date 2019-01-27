#
# This file plots the brain SPM after the GLM analysis.
#
import numpy as np
import pickle
from nilearn import plotting
import nibabel as nib

stat_list_name = "nipy_glm_stat_list.p" # Get the pickle file name of stat list.
mat_shape = (64, 64, 36) # Get the matrix shape in tuple format.
save_brain_img_name = "nipy_glm_stat_map.nii"

stat_list = pickle.load(open(stat_list_name, "rb"))
mat = np.zeros(mat_shape)

for l in stat_list:
  idx = l[1]
  val = l[0][0]
  if val > 0:
  mat[idx[0], idx[1], idx[2]] = val

brain_img = nib.Nifti1Image(mat, np.eye(4))
brain_img = nib.save(brain_img, save_brain_img_name)
brain_img = nib.load(save_brain_img_name)
plotting.plot_glass_brain(brain_img)
plotting.show()
