import os
import freesasa
import pandas as pd
import numpy as np

rootdir = "/cluster/home/reberv/picotti_nas/Viviane/Computational/AF3/all_top_models"

for file in os.listdir(rootdir):
	filename = os.fsdecode(file)
	if filename.endswith(".pdb"): 
		file_with_path = rootdir + "/" + filename
		filename_no_extension = os.path.splitext(os.path.basename(file))[0]
		filename_no_extension_with_path = rootdir + "/" + filename_no_extension
		structure = freesasa.Structure(file_with_path)
		result = freesasa.calc(structure)

		res_area = result.residueAreas()

		for key in res_area:
   		 	chain_index = np.repeat(key, len(res_area[key]))
    			aa_index = np.arange(len(res_area[key])) + 1
    
    			mainChain = []
    			sideChain = []
    			aa = []
    
    			for res_key in res_area[key]:
        			aa.append(res_area[key][res_key].residueType)
        			mainChain.append(res_area[key][res_key].mainChain)
        			sideChain.append(res_area[key][res_key].sideChain)
        
    			chain_result_dict = {'chain_index': chain_index, 'aa_index': aa_index, 'aa': aa, 'main_chain_surface': mainChain, 'side_chain_surface': sideChain}

    			chain_result_df = pd.DataFrame(data = chain_result_dict)
    
    			if key == 'A':
        			result_df = chain_result_df
    			else:
        			result_df = pd.concat([result_df, chain_result_df])
    
		result_df.to_csv(filename_no_extension_with_path + ".csv")

		continue
	else:
		continue
