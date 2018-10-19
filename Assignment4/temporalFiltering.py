#
# This file implements temporal filtering of voxel time series data.
#
# Insights into the "rest.nii.gz" data header:
#   The "dim" states: 72  72  38 120 (x, y, z, t)
#     i.e. 120 volumes (or observations) of brain is taken.
#
#   The "pixdim" states: 3.          3.          3.33000016  2.5 (x, y, z, t)
#   The "xyzt_units" states: 10
#     i.e. Each volume is taken 2.5 seconds apart i.e. time difference between
#     two samples is 2.5 seconds (or TR = 2.5 secs, also called as Time
#     Resolution parameter or sampling time).
#     In other words, the time difference between two acquisition of same slice
#     is 2.5 seconds.
#

import nibabel as nib
import numpy as np
from scipy.signal import detrend
import sys

class TemporalFiltering(object):
  def __init__(self, TR, data_path, cutoff, out_file_name, low_high=""):
    """
    Args:
      TR (float): Repetition Time in seconds (or Time Resolution).
      data_path (string): Path to the NIFTI data (including the name of file).
      cutoff (float): Cutoff time in secs i.e. Cutoff Frequency = 1.0 / cutoff.
      out_file_name (string): Output file name.
      low_high (string): <"low"|"high">
    """
    self._tr = TR
    self._data = nib.load(data_path).get_fdata()
    self._cutoff_frequency = 1.0 / cutoff
    self._low_high = low_high

    # Sampling Rate = 1.0 / TR
    self._sample_rate = 1.0 / TR # Or sample frequency.
    # Frequency Resolution of spectrum = sampling rate / number of samples
    # Number of samples = 4th (last) dimension of NIFTI data.
    self._frequency_resolution = self._sample_rate / float(self._data.shape[-1])
    self._index_of_cutoff_freq = int(
        self._cutoff_frequency / self._frequency_resolution) # Always the floor.

    # Adjust the index of cutoff frequency based on high or low pass filter.
    if self._low_high == "low": # Allow lower frequencies to pass through.
      pass
    elif self._low_high == "high": # Allow higher frequencies to pass through.
      self._index_of_cutoff_freq -= 1

    """
    TR = 2.5
    samplert = 1.0 / TR
    # Number of samples = 120
    FR = samplert / 120 # Frequency Resolution, here it is 0.0033Hz.
    """

  def _get_filtered_signal_of_voxel(self, voxel_ts):
    """
    Args:
      voxel_ts (numpy.ndarray): The timeseries array of the voxel.

    Returns:
      numpy.ndarray: A filtered timeseries signal.
    """
    voxel_ts_dt = detrend(voxel_ts)
    voxel_ts_dt_fft = np.fft.fft(voxel_ts_dt)

    if self._low_high == "high":
      # Zero down smaller or equal frequencies.
      for i in range(self._index_of_cutoff_freq):
        voxel_ts_dt_fft[i] = 0
    elif self._low_high == "low":
      # Zero down greater frequencies.
      for i in range(self._index_of_cutoff_freq, self._data.shape[-1]):
        voxel_ts_dt_fft[i] = 0

    voxel_ts_dt_fft_ifft = np.fft.ifft(voxel_ts_dt_fft)
    #return np.real(2 * voxel_ts_dt_fft_ifft)
    return np.real(voxel_ts_dt_fft_ifft)

  def do_temporal_filtering(self):
    """
    Returns:
      numpy.ndarray: A 4 Dimensional matrix.
    """
    result = np.zeros(self._data.shape)
    for x in range(self._data.shape[0]):
      for y in range(self._data.shape[1]):
        for z in range(self._data.shape[2]):
          voxel_ts = self._data[x,y,z,:]
          voxel_ts_filtered = self._get_filtered_signal_of_voxel(voxel_ts)
          result[x,y,z,:] = voxel_ts_filtered

    return result

if __name__ == "__main__":
  data_path = sys.argv[1]
  TR = float(sys.argv[2]) # in seconds.
  cutoff = float(sys.argv[3]) # in seconds.
  out_file_name = sys.argv[4]

  tempFilter = TemporalFiltering(TR, data_path, cutoff, out_file_name, "high")
  result = tempFilter.do_temporal_filtering()
  result = nib.Nifti1Image(result, np.eye(4))
  nib.save(result, ("/Personal/projects/COL786/data/data_4/sampleResult/" +
           out_file_name + "_temporal_filtered.nii.gz"))
