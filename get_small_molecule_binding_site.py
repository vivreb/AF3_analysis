import pymol2
import os


rootdir = "/cluster/home/reberv/picotti_nas/Viviane/Computational/AF3/all_top_models"

for file in os.listdir(rootdir):
	filename = os.fsdecode(file)
	if filename.endswith("atp_model.pdb"): 
		file_with_path = rootdir + "/" + filename
		export_file_with_path = rootdir + "/binding_sites/" + filename
		with pymol2.PyMOL() as pymol:
			pymol.cmd.load(file_with_path,'myprotein')
			pymol.cmd.select("binding_residues", "br. all within 5 of resn ATP")
			pymol.cmd.save(export_file_with_path.replace('atp_model.pdb', 'atp_binding_site.pdb'), selection='binding_residues')
	if filename.endswith("stu_model.pdb"): 
		file_with_path = rootdir + "/" + filename
		with pymol2.PyMOL() as pymol:
			pymol.cmd.load(file_with_path,'myprotein')
			pymol.cmd.select("binding_residues", "br. all within 5 of resn STU")
			pymol.cmd.save(export_file_with_path.replace('stu_model.pdb', 'stu_binding_site.pdb'), selection='binding_residues')