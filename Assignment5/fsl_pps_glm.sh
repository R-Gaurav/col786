#!/bin/bash
# This script is the solution for Assignment 5 of COL786. One needs to copy a
# standard design.fsf and modify it to suit ones needs. Ideally, create a
# standard design.fsf file from GUI for the first run and then copy it for next
# runs. Modification to the copied design.fsf can be automated through a bash
# script similar to this.
#
# For more information:
#   https://wiki.biac.duke.edu/biac:fsl:guide
#   http://fsl.fmrib.ox.ac.uk/fslcourse/lectures/scripting/s_0020.htm

echo "Script for FSL Preprocessing and GLM Analysis"
echo "============================================="

# Set the FMRI data path.
FD_PATH="\/Personal\/projects\/COL786\/data\/data_5\/subdata"
ANAT_IN="$FD_PATH\/anat\/sub-MSC01_ses-struct01_run-01_T1w_brain_extracted.nii.gz"
BOLD_IN="$FD_PATH\/func\/sub-MSC01_ses-func01_task-motor_run-01_bold.nii.gz"
OUTPUT="\/Personal\/projects\/COL786\/Assignment5\/output"
STCF="\/Personal\/projects\/COL786\/data\/data_5\/slice_time.txt"
IMG_REG="\/usr\/local\/fsl\/data\/standard\/MNI152_T1_2mm_brain"
declare -i N_RUNS=3
FWHM=6.0 # in mm.
ZTH=3.1
declare -i N_EVs=5

# Set the motion covariates file path.
COVR="$FD_PATH\/covariates"
EVP1="$COVR\/left_f"
EVP2="$COVR\/left_h"
EVP3="$COVR\/right_f"
EVP4="$COVR\/right_h"
EVP5="$COVR\/tongue"

# Set the explanatory variables names.
EVT1="left_f"
EVT2="left_h"
EVT3="right_f"
EVT4="right_h"
EVT5="tongue"


# Preprocessing of single subject single run data.
echo "Preprocessing single-subject single-run fMRI data..."

# The given anatomical image is already brain extracted.
echo "Anatomical Image is already brain extracted, hence no BET processing"

# View the anatomical image.
echo "Viewing the anatomical data. Close the app to proceed, when done."
#fsleyes $ANAT -cm hot

echo "Viewing the functional data. Close the app to proceed, when done."
#fsleyes $BOLD -cm hot

echo "Setting output directory @OUTPUT where FEAT output would be stored..."
sed -i '' s/@OUTPUT/$OUTPUT/g design.fsf

echo "Setting the slice time file @STCF..."
sed -i '' s/@STCF/$STCF/g design.fsf

echo "Setting the spatial smoothing @FWHM value..."
sed -i '' s/@FWHM/$FWHM/g design.fsf

echo "Setting Z threshold @ZTH..."
sed -i '' s/@ZTH/$ZTH/g design.fsf

echo "Setting up standard image registration @IMG_REG..."
sed -i '' s/@IMG_REG/$IMG_REG/g design.fsf

echo "Setting up 4D BOLD input data @BOLD_IN..."
sed -i '' s/@BOLD_IN/$BOLD_IN/g design.fsf

echo "Setting up 3D ANAT input data @ANAT_IN..."
sed -i '' s/@ANAT_IN/$ANAT_IN/g design.fsf

for i in $(seq 1 $N_EVs)
do
  echo "Setting up explanatory variable $i covariates (left_h, left_f...)..."
  evt=EVT$i
  evp=EVP$i
  sed -i '' s/@$evt/${!evt}/g design.fsf
  sed -i '' s/@$evp/${!evp}/g design.fsf
done


echo "Feat in progress..."
feat design.fsf

echo "Done"
