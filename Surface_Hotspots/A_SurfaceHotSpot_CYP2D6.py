import re
import numpy as np

####################################################  README  ###########################################################################
""" 
Purpose: Calculates cumulative number of ligand heavy atoms in a 5 Å radius of the CB atom (CA atom for glycine) of the protein amino acids. 
	 A second script ("B_SurfaceFigure_CYP2D6.py") plots the data on the surface of CYP2D6 according to a red-scale.			
																																		
Usage: Run <python3 A_SurfaceHotSpot_CYP2D6.py> in the same directory as your numbered MD frames. Make sure to adapt "user defined variables" below.		
																																		
Input: MD frames numbered from "1.pdb" to e.g. "500.pdb"							

Example: Example with 5 MD frames is provided



Written in gedit by André Fischer (2018/2019)	
"""

###############################################  USER DEFINED VARIABLES  ################################################################

ligand_type = 'ACE'				# Ligand PDB Residue Type
nfiles = 5					# Number of files for calcuation
amino_acid_atom = 'CB' 				# Atom for which radius is defined (e.g. CA, CB) - from PDB files
glycine_amino_acid_atom = 'CA'			# Atom for which radius if defined for Glycines (since no CB present)
ligand_PDB_type = 'ATOM'			# Ligand PDB Atom type
ligand_heavy_atoms = 'CNO'			# Ligand heavy atoms to be considered in alphabetic order
radius = 5					# Radius of interest
jobname = 'hotspot_example1'			# Your jobname
n_prot_res = 497				# Number of protein residues

################################################  EXECUTION  #############################################################################

#Output files
raw_filename = jobname + '.out'
raw_output = open(raw_filename, 'w')


# File loop
for i in range(1,(nfiles+1)):
		
	# Info
	print('processing file ' + str(i))
	
	# Opening file
	pdb_filename = str(i) + '.pdb'
	pdbfile_open = open(pdb_filename, 'r')
	pdbfile = pdbfile_open.read()
		
	# Header for outfile
	raw_output.write('-------File ' + str(i) + '--------')
	raw_output.write('\n')

	# Pattern for all CB atoms in PDB File
	CB_pattern = r'ATOM\s+\d+\s+' + re.escape(amino_acid_atom) + r'\s+\w+\s+\w+\s+(\d+)\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
	CB_pattern_compile = re.compile(CB_pattern)
	CB_matches = CB_pattern_compile.finditer(pdbfile)

	# Loop for CB atoms (except GLY)
	for CB_matchcount in CB_matches:
		
		# Count for residue
		per_residue_count = 0
		
		# Pattern for all ligand heavy atoms
		ligand_pattern = r'' + re.escape(ligand_PDB_type) + r'\s+\d+\s+[' + re.escape(ligand_heavy_atoms) + r']\d+\s+' + re.escape(ligand_type) + r'\s+\d+\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})' 
		ligand_pattern_compile = re.compile(ligand_pattern)
		ligand_matches = ligand_pattern_compile.finditer(pdbfile)
		
		# Loop for Ligand heavy atoms and comparison
		for ligandcounter in ligand_matches:			
			x_dist = ( float(ligandcounter.group(1)) )   -  ( float(CB_matchcount.group(2)) )
			y_dist = ( float(ligandcounter.group(2)) )   -  ( float(CB_matchcount.group(3)) )
			z_dist = ( float(ligandcounter.group(3)) )   -  ( float(CB_matchcount.group(4)) )

			squared_added = (x_dist * x_dist) + (y_dist * y_dist) + (z_dist * z_dist)
			vector_length = np.sqrt(squared_added)
			
			# Check
			if vector_length <= 5:
				per_residue_count = per_residue_count + 1
			else:
				pass	
					
		# Writing output
		raw_output.write('Residue: ' + str(CB_matchcount.group(1)) + ' | hits: ' + str(per_residue_count))
		raw_output.write('\n')
			

######### Processing Gly CA (copy of above loop)
	
	# Pattern for all GLY CA atoms in PDB File
	CA_pattern = r'ATOM\s+\d+\s+' + re.escape(glycine_amino_acid_atom) + r'\s+GLY\s+\w+\s+(\d+)\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})'
	CA_pattern_compile = re.compile(CA_pattern)
	CA_matches = CA_pattern_compile.finditer(pdbfile)	
	
	# Loop for GLY CA Atoms
	for CA_matchcount in CA_matches:
		
		# Count for residue
		per_residue_count = 0
		
		# Pattern for all ligand heavy atoms
		ligand_pattern2 = r'' + re.escape(ligand_PDB_type) + r'\s+\d+\s+[' + re.escape(ligand_heavy_atoms) + r']\d+\s+' + re.escape(ligand_type) + r'\s+\d+\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})\s+([0-9-]+\.\d{3})' 
		ligand_pattern_compile2 = re.compile(ligand_pattern2)
		ligand_matches2 = ligand_pattern_compile2.finditer(pdbfile)
		
		# Loop for Ligand heavy atoms and comparison
		for ligandcounter2 in ligand_matches2:		
			x_dist = ( float(ligandcounter2.group(1)) )   -  ( float(CA_matchcount.group(2)) )
			y_dist = ( float(ligandcounter2.group(2)) )   -  ( float(CA_matchcount.group(3)) )
			z_dist = ( float(ligandcounter2.group(3)) )   -  ( float(CA_matchcount.group(4)) )

			squared_added = (x_dist * x_dist) + (y_dist * y_dist) + (z_dist * z_dist)
			vector_length = np.sqrt(squared_added)
			
			# Check
			if vector_length <= 5:
				per_residue_count = per_residue_count + 1
				
			else:
				pass
						
		# Handling output
		raw_output.write('Residue: ' + str(CA_matchcount.group(1)) + ' | hits: ' + str(per_residue_count))
		raw_output.write('\n')

	# Closing raw output file
raw_output.close()

######### Output post-processing
print('--Processing output--')
procme_open = open(raw_filename, 'r')
procme_file = procme_open.read()

valuefile_name = jobname + '.res'
valuefile = open(valuefile_name, 'w')
valuefile.write('Calculation: ' + jobname)
valuefile.write('\n')

resultfile_name = jobname + '_clean.res'
resultfile = open(resultfile_name, 'w')
resultfile.write('Calculation: ' + jobname)
resultfile.write('\n')


# Residue loop
for resnr in range(1,(n_prot_res + 1)):
	current_residue_pattern = r'Residue\:\s+' + re.escape(str(resnr)) + r'\s+\|\s+hits\:\s+(\d+)'
	current_residue_pattern_compile = re.compile(current_residue_pattern)
	current_residue_matches = current_residue_pattern_compile.finditer(procme_file)
	
	residue_hits_count = 0
	
	for resnr_c in current_residue_matches:
		residue_hits_count = residue_hits_count + int(resnr_c.group(1))
		
	
	# Output for residue
	valuefile.write('Residue ' + str(resnr) + ' | ' + str(residue_hits_count))
	valuefile.write('\n')
	
	# Clean Output (only residues with hits)
	if residue_hits_count > 0:
		print('Residue ' + str(resnr) + ' | ' + str(residue_hits_count))
		resultfile.write('Residue ' + str(resnr) + ' | ' + str(residue_hits_count))
		resultfile.write('\n')
	else:
		pass

