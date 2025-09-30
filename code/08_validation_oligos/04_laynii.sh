#!/bin/bash

d0=$PWD

FSLOUTPUTTYPE=NIFTI

#for sub in Br5766; do echo "+++ $sub"
for sub in Br2934 Br5347 Br5394 Br5766 Br5785 Br9413; do echo "+++ $sub"

	subImgDir=$d0/../../plots/08_validation_oligos/$sub
	
	echo "make layers"
	LN2_LAYERS -rim $subImgDir/${sub}_128x128binned_manual_compartments.nii.gz \
		-nr_layers 4 -equivol -iter_smooth 5000 \
		-output $subImgDir/${sub}_128x128binned_laynii

	echo "define SWM mask"
	# 0.4792um is 1x1 bin resolution
	# layers file has resolution listed as 1mm x 1mm
	# therefore, to define 1mm SWM
	# 1mm * 1000um/(0.4792*128) = 15.71
	LN2_CHOLMO -layers $subImgDir/${sub}_128x128binned_laynii_layers_equivol.nii \
		-outer -nr_layers 1 -layer_thickness 15.71 \
		-output $subImgDir/${sub}_128x128binned_laynii_SWM.nii
	fslmaths $subImgDir/${sub}_128x128binned_laynii_SWM.nii \
		-thr 5 -bin \
		$subImgDir/${sub}_128x128binned_laynii_SWM.nii
	

	for ii in equivol equidist; do echo "+++ $ii"
		
		echo "intersect with tissue mask"
		fslmaths $subImgDir/${sub}_128x128binned_laynii_layers_${ii}.nii \
			-mul $subImgDir/${sub}_128x128binned_manual_tissue_mask.nii \
			$subImgDir/${sub}_128x128binned_laynii_layers_${ii}_masked.nii

	done
	
	echo "intersect white matter with tissue mask"
	fslmaths $subImgDir/${sub}_128x128binned_laynii_SWM.nii \
		-mul $subImgDir/${sub}_128x128binned_manual_tissue_mask.nii \
                $subImgDir/${sub}_128x128binned_laynii_SWM_masked.nii

done
