#!/bin/bash

# script to facilitate manual drawing of GM/WM & pial tissue boundaries

d0=$PWD
sub=$1
subImgDir=$d0/../../plots/08_validation_oligos/$sub

FSLOUTPUTTYPE=NIFTI

echo ""
echo "Save ${sub}_128x128binned_manual_boundaries with:"
echo "GM/WM 	= 1"
echo "pial	= 2"
echo ""
echo "Save ${sub}_128x128binned_manual_compartments with:"
echo "WM 	= 1"
echo "GM 	= 3"
echo "CSF 	= 2"
echo ""
echo echo "Save ${sub}_128x128binned_manual_white_matter with:"
echo "WM    = 1"
echo ""
echo echo "Save ${sub}_128x128binned_manual_tissue_mask with:"
echo "tissue    = 1"
echo ""

fslview_deprecated -m single \
	$subImgDir/${sub}_array_DAPI_DAPI_Positive_128x128binned_for_masking.nii -l "Red" -b 0,5 \
	$subImgDir/${sub}_array_MBP_Oligo_128x128binned_for_masking.nii -l "Blue-Lightblue" -b 0,5 \
	$subImgDir/${sub}_array_HPCAL1_L2_neuron_128x128binned_for_masking.nii -l "Red-Yellow" -b 0,5 
