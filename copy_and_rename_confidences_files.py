import os

rootdir = "/cluster/scratch/reberv/json_files_kinases"

for f in os.scandir(rootdir):
	if f.is_dir():
		subdir1 = f.path + "/alphafold3_predictions"
		if(os.path.exists(subdir1)):
			for f2 in os.scandir(subdir1):
				if f2.is_dir():
					subdir2 = f2.path
					for file in os.listdir(subdir2):
						filename = os.fsdecode(file)
						if filename.endswith("confidences.json"): 
							if filename.endswith("summary_confidences.json"): 
								continue
							else:
								cmd2 = "cp " + subdir2 + "/" + filename + " /cluster/scratch/reberv/all_confidences/"
								os.system(cmd2)
								continue