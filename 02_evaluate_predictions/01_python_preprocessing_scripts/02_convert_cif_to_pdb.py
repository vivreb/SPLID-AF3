import pymol2
import os


rootdir = "/cluster/scratch/reberv/all_top_cif_files"

for file in os.listdir(rootdir):
	filename = os.fsdecode(file)
	if filename.endswith(".cif"): 
		file_with_path = rootdir + "/" + filename
		with pymol2.PyMOL() as pymol:
     			pymol.cmd.load(file_with_path,'myprotein')
     			pymol.cmd.save(file_with_path.replace('.cif', '.pdb'), selection='myprotein')