import matplotlib.cm as cm
import matplotlib.colors as colors
from pymol import cmd
import re

####################################################  README  ###########################################################################

""" 
Purpose: Script is designed to be used after hotspot calculation using "A_SurfaceHotSpot_CYP2D6.py" is finished. It plots the hotspots on
the surface of CYP2D6 protein that is provided within Pymol Session file.
																																		
Usage: Open Pymol, load session "Surface_Pymol_Session.pse". Within Python, run the python script "B_SurfaceFigure_CYP2D6.py" to display hotspots on the surface		
																																		
Input: Output file from hotspot calculation (not "jobname_clean.res", but "jobname.res")



Written in gedit by André Fischer (2018/2019)	
"""

###############################################  USER DEFINED VARIABLES  ################################################################

nres = 497						# Number of protein residues
input_file = 'hotspot_example1.res'			# Name of input file from hotspot calculation
color_map = 'Reds'  					# Others: 'cool_warm',...
lower_value_range = 0					# Lowest value that should get colored
upper_value_range = 5					# Highest value that should get colored

################################################  EXECUTION  #############################################################################

norm = colors.Normalize(vmin=lower_value_range, vmax=upper_value_range)
f2rgb = cm.ScalarMappable(norm=norm, cmap=cm.get_cmap(color_map))
def f2hex(f2rgb, f):
    rgb = (f2rgb.to_rgba(f)[:3])
    return rgb

# Opening input file
resfile_open = open(input_file, 'r')
resfile = resfile_open.read()
	
# Residue loop
for i in range(1,(nres+1)):
	current_residues_pattern = r'Residue\s+' + re.escape(str(i)) + r'\s+\|\s+(\d+)'
	crp_compile = re.compile(current_residues_pattern)
	current_matches = crp_compile.finditer(resfile)
	
	# Looping through matches
	for match_counter in current_matches:
		x = int(match_counter.group(1))
	
	print(x)	# Number of hits
	
	v1 = f2hex(f2rgb, x)
	v1_r = v1[0] * 255
	v1_g = v1[1] * 255
	v1_b = v1[2] * 255
	
	color_tupple = (v1_r,v1_g,v1_b)
	
	print('RGB Tuppel for residue ' + str(i))
	print(v1_r,v1_g,v1_b)


	# Creating color for pymol
	colname = 'col' + str(i)
	current_color = cmd.set_color(colname, color_tupple)
	
	# Selecting corrrect amino acid
	resi_string = 'resi ' + str(i)
	resi_select = 'sel' + str(i)
	
	# Selection of correct AA
	cmd.select(resi_select, resi_string)
	
	# Coloring of correct AA
	cmd.color(colname, resi_select)
	






