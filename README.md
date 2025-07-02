# Analysis of structural changes as a result of small-molecule inhibitor binding

We use a split-domain approach to generate structural predictions with and without a small-molecule ligand. The resulting structural changes are analysed using the presented code.

## How to use

Save "generate_AF3_input_json_files.R" and the file "231129_human_kinome.xlsx" into the same folder. In this folder, create a subdirectory called "json_files_kinases". The output from "generate_AF3_input_json_files.R" will be saved in this file. Run "generate_AF3_input_json_files.R" line by line. The output might need to be modified if the AlphaFold3 predictions are not according to https://github.com/jurgjn/batch-infer on the ETH Euler cluster.

Next, jobs are submitted by running "setup_jobs_script.py" on the cluster. It takes approximately 24 hours for all predictions to run.

Results are copied by running "copy_and_rename_pdb_files.py" and "copy_and_rename_confidences_files.py" on the cluster. The directories containing these files can be copied onto another server at this point.

Next, the exported .cif files are converted to .pdb files by running "convert_cif_to_pdb.py" on the cluster. From these .pdb files, predicted structural features can be extracted. This is done by running "extract_confident_residues.py", "get_small_molecule_binding_site.py" and "get_surface_area_with_freesasa.py", to extract regions of high confidence, regions that bind the ligand, and regions that are on the surface, respectively.

Make sure that in the end, you have the following folder structure: in your working directory, the folder "all_summary_confidences" contains all generated "[protein_id]_summary_confidences.json" files, and the folder "all_top_models" contains the .cif and .pdb files. Within the "all top models" folder, you should also have the surface accessibility saved in .csv files, and the two folders "binding_sites" and "confident_residues", containing the .pdb files of the predicted small molecule binding sites and the high confidence residues, respectively. 

In this working directory, save the scripts "functions_for_AF3_structural_analysis.R", "AF3_structural_analysis.R", and "set_cutoffs_for_RMSSD_and_SASA.R. Run "AF3_structural_analysis.R" line by line, followed by "set_cutoffs_for_RMSSD_and_SASA.R" line by line. 

The output is a file containing a list of all the proteins with all the quantities calculated in the analysis.