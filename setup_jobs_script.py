import os

rootdir = "/cluster/scratch/reberv/json_files_kinases"

for f in os.scandir(rootdir):
	if f.is_dir():
		subdir = f.path
		cmd2 = "./batch-infer alphafold3 " + subdir + " | sbatch"
		os.system(cmd2)