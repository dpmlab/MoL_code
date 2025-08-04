#first variable is BIDS subject label

	# -B /data/MoL_clean/BIDS/:/data \

singularity run --cleanenv \
	-B /data/MoL_experts/data/:/data \
	-B /data/MoL_clean/fmriprep:/out \
	-B /data/MoL_clean/scratch:/work \
	fmriprep-23.0.2.simg \
	/data /out \
	participant \
	--ignore=slicetiming \
	--use-syn-sdc \
	--fs-license-file=/out/fs-license.txt \
	--output-spaces MNI152NLin2009cAsym fsaverage6 \
	-w /work/ \
	--low-mem \
	--participant-label $1 

