#!/bin/bash

d0=$PWD

FSLOUTPUTTYPE=NIFTI

for sub in Br5766; do echo "+++ $sub"
#for sub in Br2934 Br5347 Br5394 Br5766 Br5785 Br9413; do echo "+++ $sub"

	subImgDir=$d0/../../plots/08_validation_oligos/$sub
	
	# make layers
	LN2_LAYERS -rim $subImgDir/${sub}_128x128binned_manual_compartments.nii.gz \
		-nr_layers 4 -equivol -iter_smooth 5000 \
		-output $subImgDir/${sub}_128x128binned_laynii

	for ii in equivol equidist; do echo "+++ $ii"
		
		echo "intersect with tissue mask"
		fslmaths $subImgDir/${sub}_128x128binned_laynii_layers_${ii}.nii \
			-mul $subImgDir/${sub}_128x128binned_manual_tissue_mask.nii \
			$subImgDir/${sub}_128x128binned_laynii_layers_${ii}_masked.nii

	done
	
	echo "intersect white matter with tissue mask"
	fslmaths $subImgDir/${sub}_128x128binned_manual_white_matter.nii \
		-mul $subImgDir/${sub}_128x128binned_manual_tissue_mask.nii \
                $subImgDir/${sub}_128x128binned_manual_white_matter_masked.nii

done
