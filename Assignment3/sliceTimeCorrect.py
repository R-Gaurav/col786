#
# Slice Time Correction Code.
#
# There are three popular methods to do Slice Time Correction (STC), namely:
# linear sampling, cubic spline and sinc sampling. This file implements linear
# sampling method. However linear sampling results in smoothing of image, which
# varies as per the time point at which sampling is done. That is: For a time
# point (at which the data is to be sampled) closer to either of the actual
# measured time points, smoothing is less (and this is preferred also the result
# can be expected to that obtained from either cubic spline or sinc methods),
# but if the sampling point is in between the measured time points then smoothing
# is rather strong and this leads to loss of information, hence not preferred.
#
# However with the case of cubic spline or sinc methods, choice of sampling time
# point is not much of issue, the results are same.
#
# In this data of COL786, 8 slices are excited are excited simultaneously and
# are acquired in 9 batches with interleaved acquisition time. In other words
# the acceleration factor is 8.
#
# Links:
# http://brainvoyager.com/bv/doc/UsersGuide/Preprocessing/SliceScanTimeCorrection.html
# https://matthew-brett.github.io/teaching/slice_timing.html
#

import nibabel as nib
import numpy as np
import sys

def _read_slice_time_acquisition_file_list(file_path):
  """
  Reads the file in passed `file_path` and returns a list of floating point
  numbers which denote the time at which sliced were acquired on Z axis.

  Note: The length of returned list should be same as the number of slices along
  Z axis. The float at index 0 is the time at which slice 0 (at the bottom of
  Z axis was acquired) and so on... i.e. slices are numbered from bottom of Z
  axis.

  Args:
    file_path (str): Path to the slice time acquisition file.

  Returns:
    [float]
  """
  f = open(file_path, "r")
  lines = f.readlines()
  f.close()
  slice_times = [float(st.rstrip("\r\n")) for st in lines]
  return slice_times

def _is_slice_timings_valid(slice_times, TR):
  """
  Validates slice times in `slice_times` list.

  Args:
    slice_times (list): A list of floating points which indicate the time when
        each slice is captured.
    TR (float): Repeatition Time.

  Returns:
    bool
  """
  for st in slice_times:
    if st > TR: # If any slice time is greater than TR, then the list is invalid.
      return False

  return True

def _write_file(sof, msg, output_file):
  """
  Writes `sof` to `file_name`.txt

  Args:
    sof (str): <"SUCCESS"|"FAILURE">
    msg (str): Message to be printed.
    output_file (str): Output file name.
  """
  f = open(output_file+".txt", "w")
  f.write(sof+"\n")
  print msg
  f.close()
  sys.exit()

def _2d_linear_interpolate(fst_vol_slice, lst_vol_slice, fv_slc_tme, lv_slc_tme,
    target_time):
  """
  It linearly interpolates from the first and last vol slice at the target time.

  Args:
    fst_vol_slice (numpy.ndarray): A 2D matrix (Y1).
    lst_vol_slice (numpy.ndarray): A 2D matrix (Y2).
    fv_slc_tme (float): Time when first slice was taken (X1).
    lv_slc_tme (float): Time when last slice was taken (X2).
    target_time (float): A target time lying between the times of first slice
        and last slice (X).
  """
  # Find the slope of linear line from (X1, Y1) to (X2, Y2).
  m = (lst_vol_slice - fst_vol_slice) / (lv_slc_tme - fv_slc_tme)
  return m * (target_time - fv_slc_tme) + fst_vol_slice

def do_slice_time_correction(img_file, TR, target_time, st_acq_file,
    output_file):
  """
  Does slice time correction of the given NIFTI img_file. Saves the output to a
  a new file with name `output_file`.nii.gz. It also outputs "SUCCESS" and
  "FAILURE" based on the input.

  Note: If the target_time is between lower bound and upper bound of image
  acquisition time then it is SUCCESS else FAILURE. Also all the slice
  acquisition time in `st_acq_file` should be less than or equal to TR.

  Args:
    img_file (str): Path/to/the/NIFTI/*.nii.gz.
    TR (float): Repeatition Time, i.e. the time between two consecutive RF pulse.
    target_time (float): The time to which all the slices are to be corrected,
        i.e. it should look like that the slices in the functional volume
        (brain) were all acquired simulatneously at the target time.
    st_acq_file (str): Path/to/the/sliceTimeAcquisitionFile.txt
    output_file (str): Output file name.

  Returns:
    numpy.ndarray(:,:,:,:): A 4D array.
  """
  slice_times = _read_slice_time_acquisition_file_list(st_acq_file)

  # Find if all the slice times are valid or not.
  if not _is_slice_timings_valid(slice_times, TR):
    _write_file("FAILURE", "Slice times are incorrect.", output_file)

  # Find if the `target_time` is valid or not.
  lower_bound = min(slice_times)
  upper_bound = max(slice_times)
  if not (lower_bound <= target_time and target_time <= upper_bound):
    _write_file("FAILURE", "Target Time is out of [%s, %s] bound" % (
                lower_bound, upper_bound), output_file)

  img = nib.load(img_file)
  data = img.get_fdata()
  shape = data.shape # This returns a 4D tuple: (90, 104, 72, 100).

  # The data in the assignment has 72 slices and 100 brain volumes i.e. the
  # dimension of the time series is 100.

  # Check if the length of `slice_times` list is equal to the Z dimension.
  if not (len(slice_times) == shape[2]):
    _write_file("FAILURE", ("Length of slice times (: %s) is not equal to the Z"
                "dimension (: %s)" % (len(slice_times), shape[2])), output_file)

  # Proceed to implement the linear sampling at the `target_time`.
  # Note that first and last volume i.e. data[:,:,:,0] and data[:,:,:,shape[3]-1]
  # are supposed to be constant as mentioned in assignment, hence no special
  # trick for linear interpolation is required for them.
  # Every slice in each volume (remaining volumes after removing first and last
  # ones) needs to be linearly interpolated.
  # For the sake of convenience, Z axis slices are numbered from [0 to 71].

  result = np.ndarray(shape)
  result[:,:,:, 0] = data[:,:,:, 0]
  # For every volume.
  for crnt_vol in range(1,shape[3]-1):
    prev_vol = crnt_vol - 1
    next_vol = crnt_vol + 1

    # For every slice in each current volume do linear interpolation.
    for slc_num in range(shape[2]):
      crnt_vol_slc_time = slice_times[slc_num]
      crnt_vol_slc = data[:,:, slc_num, crnt_vol]

      # Depending on the target time and current volume slice time, we need to
      # figure out that whether we need a previous volume slice or next volume
      # slice to the current one for linear interpolation.
      if crnt_vol_slc_time == target_time:
        # There is no need to linearly interpolate this slice as this is already
        # at the desired target time.
        result_slice = crnt_vol_slc

      if crnt_vol_slc_time < target_time:
        next_vol_slc = data[:,:, slc_num, next_vol] # Next volume slice.
        result_slice = _2d_linear_interpolate(
            crnt_vol_slc, next_vol_slc, crnt_vol_slc_time,
            crnt_vol_slc_time + TR, target_time)
      if crnt_vol_slc_time > target_time:
        prev_vol_slc = data[:,:, slc_num, prev_vol] # Previous volume slice.
        result_slice = _2d_linear_interpolate(
            prev_vol_slc, crnt_vol_slc, crnt_vol_slc_time - TR,
            crnt_vol_slc_time, target_time)

      result[:,:, slc_num, crnt_vol] = result_slice

  result[:,:,:, shape[3]-1] = data[:,:,:, shape[3]-1]
  return result

if __name__ == "__main__":
  img_file = sys.argv[1]
  TR = float(sys.argv[2])
  target_time = float(sys.argv[3])
  st_acq_file = sys.argv[4]
  output_file = sys.argv[5]
  result = do_slice_time_correction(
      img_file, TR, target_time, st_acq_file, output_file)
  # Provide image data and an image coordinate transformation (affine).
  img = nib.Nifti1Image(result, np.eye(4))
  nib.save(
      img, ("/Personal/projects/COL786/data/data_3/example/" + output_file +
      ".nii.gz"))
  _write_file(
      "SUCCESS", "The given volume is now slice time corrected.", output_file)
