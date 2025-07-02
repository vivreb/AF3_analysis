import pymol2
import os


rootdir = "/cluster/home/reberv/picotti_nas/Viviane/Computational/AF3/all_top_models"

for file in os.listdir(rootdir):
	filename = os.fsdecode(file)
	if filename.endswith(".pdb"): 
		file_with_path = rootdir + "/" + filename
		export_file_with_path = rootdir + "/confident_residues/" + filename
		with pymol2.PyMOL() as pymol:
			pymol.cmd.load(file_with_path,'myprotein')
			pymol.cmd.select("confident_residues", "b > 70")
			pymol.cmd.save(export_file_with_path.replace('.pdb', '_confident.pdb'), selection='confident_residues')