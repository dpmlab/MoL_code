# MoL_code
Code for https://www.biorxiv.org/content/10.1101/2024.12.19.629352v2
The Neuroimaging data can be found in doi:10.18112/openneuro.ds005894.v1.0.0


Updated sheets: folder with excel sheets with cleaned loci description and encoding and retrieval, what is spoken and their respective timing. Used for creating design matrices in GLM and labeling representations with corresponding locus/item/combined representations. Spoken recall is used for semantic analyses. It also has a sheet recording the performance at screening.

performance_scoring.ipynb calculates and visualize participants’ performance scores

Run fMRI prep on the BIDS files. (Using runSingularity.sh - need to have Singularity container for the right version of fMRIprep)

After fMRI prep, run post_fmriprep.ipynb to remove nuisance regressors, get surface data + hippocampus data (with both individual masks and group masks). The output is saved in the “preprocessed” folder. Subsequent analyses use the .h5 files generated in this process.

Use get_beta_patterns.ipynb to extract multivariate beta patterns for locus, item, encoding, and retrieval. 

pat_sim_analysis.ipynb generates loci-loci, item-item, loci-encoding, item-encoding representational similarity, and weights of locus, item, encoding residuals predicting retrieval. 

hipp_sem.ipynb creates sheets where each row is a locus-item pair, that include information on representational similarity between task runs, regression weights (in predicting retrieval), & univariate activity in different ROIs, also story deviation, whether recall is correct, speak duration, etc… Create a dataframe used for subsequent analyses in R, with analysis_notebook.rmd and rsa.ipynb

rsa.ipynb looks at w4d2 scan, showed how semantic similarity in stories is related to neural similarity (from different participants for the same locus-item pair). 

visualization_mac.ipynb plot surface maps for the pat_sim_analysis 
