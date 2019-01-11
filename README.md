This dir contains the assignment details of [COL 786: Advanced Functional
Brain Imaging course](http://www.cse.iitd.ernet.in/~rahulgarg/Teaching/2016/COL786.htm).

The course consisted of 6 assignments, namely:

* Assignment 1:
Given a T1 and a T2 weighted image, identify four major lobes of brain and the
distinguishing landmarks separating those lobes.

* Assignment 2:
Given a functional image and a voxel coordinates (x, y, z), write a program to
print the time series of the voxel. [Link to code](https://github.com/R-Gaurav/col786/blob/master/Assignment2/voxel_time_series_generator.py)

* Assignment 3:
Given a functional image, TR, a slice time acquisition file which contains the
time at which each slice was acquired and target time, write a program to
perform slice time correction using linear interpolation. [Link to code](https://github.com/R-Gaurav/col786/blob/master/Assignment3/sliceTimeCorrect.py)

  Note: The slices are acquired in parallel.
  
  Assumption: The first and last volume are kept constant.

* Assignment 4:

  * 4a) Given a functional image and FWHM (in mm), perform spatial smoothing.
  [Link to code](https://github.com/R-Gaurav/col786/blob/master/Assignment4/spatialSmoothing3D.py)

  * 4b) Given a functional image, TR and cut-off time, perform temporal
  filtering. [Link to code](https://github.com/R-Gaurav/col786/blob/master/Assignment4/temporalFiltering.py)

  Note: Only FFT and IFFT libraries can be used for above both.

* Assignment 5:
Given a task based functional data (Buckner et al. 2011; Yeo et al. 2011) of a
subject, his structural data, five covariates one for each of the five motor
tasks and slice time acquisition file, use FSL tool via command line interface
to perform all the preprocessing and GLM analysis.

  Note: The subject is presented with visual cues to either tap their left or right
  fingers, squeeze their left or right toes or move their tongue to map his motor
  areas. Each block of a movement type lasts 12 seconds (10 movements) and is
  preceded by a 3 seconds cue. In each run there are 13 blocks, with 2 of tongue
  movements, 4 of hand movements (2 right, 2 left) and 4 of foot movements (2
  right, 2 left). In addition, there are 3 15-seconds fixation blocks per run.

  Find contrasts for all five motor tasks. Use 3mm for spatial smoothing and 100s
  high pass filter. [Link to code](https://github.com/R-Gaurav/col786/tree/master/Assignment5)

* Assignment 6:
  With respect to Assignment 5, perform preprocessing using FSL tool and write a
  program for performing GLM analysis to find z-stat values. Use FSL again for
  full model processing on same data and compare the z-stat values. [Link to code](
  https://github.com/R-Gaurav/col786/blob/master/Assignment6/python_glm.py)

  Note: Perform analysis for only the left foot movement covariate.
