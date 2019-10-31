# coding=utf-8
import re
import os
import numpy as np
import matplotlib.pyplot as plt


####################################################  README  ###########################################################################
""" 
Purpose: Calculates distances between mass center of all ligands in simulation to calculate average + SD; outputs raw data + plot 		

Usage: Run <python3 Aggregation_CYP2D6.py> in the same directory as your MD frames. Make sure to adapt "user defined variables" below.

Input: MD frames numbered for "os" module to read sequentially

Example: Example with 99 MD frames is provided



Written in PyCharm by André Fischer (2019)										

"""

###############################################  USER DEFINED VARIABLES  ################################################################
residue_name = 'DEB'                        # residue name specified in PDB file
number_of_residues = 20                     # number of residues in simulation
identifier = 'DEB2'                         # identifier of the simulation (random name)
starting_number = 901                       # starting number of first ligand in PDB file
heavy_atom_types = 'CNO'                    # types of ligand heavy atoms to be considered (alphabetically); 'CNO' = carbon, nitrogen, oxygen
first_frame = 0                             # first frame (for plotting)
times_inclusion = True                      # use simulation time instead of frames for plotting
frame_timestep = 0.048                      # simulation timestep in nanoseconds (ns) (only needed if time_inclusion is True)
plot_limits = False                         # make the plot start in bottom corner if "False", if "True", there will be a gap
max_simulation_time = 4.752                 # simulation time covered by the frames
x_ticks = [0,2.376,4.752]                   # ticks for the plot (simulation time)
boundary_box = [115.369,113.916,158.816]    # boundary box size x,y,z

################################################  EXECUTION  #############################################################################

def pdbReader(path):
    """
    Looks for PDB Files in supplied directory; returns sorted array with file names
    """
    pdb_file_array = []
    for r, d, f in os.walk(path):
        for file in f:
            if '.pdb' in file:
                pdb_file_array.append(file)
    pdb_file_array.sort()

    return pdb_file_array

def pointDistance(ax,ay,az,bx,by,bz):
    """
    Returns the distance between two points with coords P1(ax,ay,az) and P2(bx,by,bz) based
    on vector geometry; returns float distance
    """
    dist_x = float(ax) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az) - float(bz)

    dist_x_sq = dist_x*dist_x
    dist_y_sq = dist_y*dist_y
    dist_z_sq = dist_z*dist_z

    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    return v_length

def pointDistanceBoundary(ax,ay,az,bx,by,bz, chorus_x, chorus_y, chorus_z):
    """
    Returns the distance between two points with coords P1(ax,ay,az) and P2(bx,by,bz) based
    on vector geometry; returns float distance, but accounts for peridoic boundary conditions

    """
    v_length_array = []

    # Regular Box
    dist_x = float(ax) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az) - float(bz)

    dist_x_sq = dist_x*dist_x
    dist_y_sq = dist_y*dist_y
    dist_z_sq = dist_z*dist_z

    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # translators:
    ax_1 = ax + chorus_x
    ax_2 = ax - chorus_x
    ay_1 = ay + chorus_y
    ay_2 = ay - chorus_y
    az_1 = az + chorus_z
    az_2 = az - chorus_z

    # (0, 0, +z)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x*dist_x
    dist_y_sq = dist_y*dist_y
    dist_z_sq = dist_z*dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (0, 0, -z)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x*dist_x
    dist_y_sq = dist_y*dist_y
    dist_z_sq = dist_z*dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+x, 0, 0)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (-x, 0, 0)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+x, 0, +z)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (-x, 0, +z)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+x, 0, -z)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (-x, 0, -z)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+x, +y, +z)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (-x, +y, +z)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+x, -y, +z)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (-x, -y, -z)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (0, +y, 0)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (0, -y, 0)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (0, +, -)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # ( 0 + +)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (0 - - )
    dist_x = float(ax) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+ + -)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+ + 0)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (- - 0)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (- + 0)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+ - 0)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (-, +, -)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay_1) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (- - +)
    dist_x = float(ax_2) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (+ - -)
    dist_x = float(ax_1) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az_2) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))

    # (0 - +)
    dist_x = float(ax) - float(bx)
    dist_y = float(ay_2) - float(by)
    dist_z = float(az_1) - float(bz)
    dist_x_sq = dist_x * dist_x
    dist_y_sq = dist_y * dist_y
    dist_z_sq = dist_z * dist_z
    v_length = np.sqrt(dist_x_sq + dist_y_sq + dist_z_sq)
    v_length_array.append(float(v_length))
    v_min_length = min(v_length_array)

    return v_min_length

#////MAIN

if __name__=='__main__':

    # array for average distance between ligands per frame (+sd)
    average_frame_distances_array = []
    stdev_frame_distances_array = []

    # reading pdbs in working directory
    working_directory = os.getcwd() + '/'
    pdb_file_array = pdbReader(working_directory)

    # creating array with residue numbers (to not use range below)
    resnr_array = []
    number_of_residues = number_of_residues + starting_number
    for nr in range(starting_number, number_of_residues):
        resnr_array.append(int(nr))
    print('Accessing residue "' + residue_name + '" '  + str(len(resnr_array))  + 'x times with the following numbers: ')
    print(resnr_array)

    # opening output file
    aggregation_output_file = open('aggreation_' + identifier + '.out', 'w')

    # array to count frames for time array generation below
    frame_counter = first_frame
    frame_counter_array = []

    print('calcualting distances...')
    # Looping through files
    for pdb_file in pdb_file_array:
        print('processing file: ' + pdb_file + '...')
        frame_counter_array.append(int(frame_counter))
        frame_counter = frame_counter + 1   # incrementing framecounter

        # Array to fill with distances for each PDB frame
        distances_in_frame_array = []

        # reading file
        current_pdb_open = open(pdb_file,'r')
        current_pdb = current_pdb_open.read()

        # Looping through residue for calculation
        for outer_resnr in resnr_array:
            ## getting mass center of current ligand
            # getting coords
            outer_x_coords = []
            outer_y_coords = []
            outer_z_coords = []

            HA_pattern = r'\S+\s+[' + re.escape(heavy_atom_types) + r']\S+\s+' + re.escape(residue_name) + r'\s+' + re.escape(str(outer_resnr)) + r'\s+(\S+)\s+(\S+)\s+(\S+)'
            HA_pattern_c = re.compile(HA_pattern)
            HA_matches = HA_pattern_c.finditer(current_pdb)
            for m1 in HA_matches:
                outer_x_coords.append(float(m1.group(1)))
                outer_y_coords.append(float(m1.group(2)))
                outer_z_coords.append(float(m1.group(3)))


                if pdb_file == pdb_file_array[0] and outer_resnr == resnr_array[0]:
                    print(m1.group(0))

            # averaging coords
            outer_x_MC = (sum(outer_x_coords)) / (len(outer_x_coords))
            outer_y_MC = (sum(outer_y_coords)) / (len(outer_y_coords))
            outer_z_MC = (sum(outer_z_coords)) / (len(outer_z_coords))
            #print('outer residue MC: ' + str(outer_x_MC) + ' ' + str(outer_y_MC) + ' ' + str(outer_z_MC))

            # Looping through all other residues
            for inner_resnr in resnr_array:
                # only execute if residue numbers are different
                if inner_resnr != outer_resnr:

                    # coord arrays
                    inner_x_coords = []
                    inner_y_coords = []
                    inner_z_coords = []

                    # calculating mass center of other ligand
                    inner_HA_pattern = r'\S+\s+[' + re.escape(heavy_atom_types) + r']\S+\s+' + re.escape(residue_name) + r'\s+' + re.escape(str(inner_resnr)) + r'\s+(\S+)\s+(\S+)\s+(\S+)'
                    inner_HA_pattern_c = re.compile(inner_HA_pattern)
                    inner_HA_matches = inner_HA_pattern_c.finditer(current_pdb)
                    for m2 in inner_HA_matches:
                        inner_x_coords.append(float(m2.group(1)))
                        inner_y_coords.append(float(m2.group(2)))
                        inner_z_coords.append(float(m2.group(3)))

                    # averaging inner residue coords to get MassCenter
                    inner_x_MC = (sum(inner_x_coords)) / (len(inner_x_coords))
                    inner_y_MC = (sum(inner_y_coords)) / (len(inner_y_coords))
                    inner_z_MC = (sum(inner_z_coords)) / (len(inner_z_coords))
                    #print('inner residue MC: ' + str(inner_x_MC) + ' ' + str(inner_y_MC) + ' ' + str(inner_z_MC))

                    # calcualting distance between the two MCs
                    current_distance = pointDistanceBoundary(outer_x_MC, outer_y_MC, outer_z_MC, inner_x_MC, inner_y_MC, inner_z_MC, boundary_box[0], boundary_box[1], boundary_box[2])

                    distances_in_frame_array.append(current_distance)

                else:
                    #print('hitting exception, not calcualting for ' + str(inner_resnr))
                    pass

        # calcuating average of all distances that were measured in current frame
        average_frame_distance = (sum(distances_in_frame_array)) / (len(distances_in_frame_array))
        average_frame_distances_array.append(float(average_frame_distance))
        aggregation_output_file.write('File ' + str(pdb_file) + ' Distance: ' + str(average_frame_distance))
        aggregation_output_file.write('\n')

        # numpy stdev
        np_distances_in_frame_array = np.array(distances_in_frame_array)
        np_stdev_average_frame = np.std(np_distances_in_frame_array)
        stdev_frame_distances_array.append(float(np_stdev_average_frame))

    #### Plotting
    print('plotting data...')
    # "Preparing" data
    plot_x = frame_counter_array
    plot_y = average_frame_distances_array
    x_label = 'Frames'

    ### Standard deviation
    plot_sd_array = stdev_frame_distances_array

    # creating new datasets representing SD-lines
    plot_y_plus_std = []
    plot_y_minus_std = []

    stdev_array_counter = 0
    for v1 in plot_y:
        value_plus = v1 + plot_sd_array[stdev_array_counter]
        value_minus = v1 - plot_sd_array[stdev_array_counter]

        stdev_array_counter = stdev_array_counter + 1

        plot_y_plus_std.append(value_plus)
        plot_y_minus_std.append(value_minus)

    # Getting array with simulation time instead of frame numbers
    if times_inclusion == True:         # only executes if user specifies it above
        times_array = []
        for frame_number in frame_counter_array:
            current_time = frame_number * frame_timestep
            times_array.append(float(current_time))
        plot_x = times_array    # replacing frame array with time array for x-coord of plot
        x_label = 'Simulation time (ns)'

    # Creating a matplotlib plot and saving as PNG
    fig, ax1 = plt.subplots()
    color = '#006e6e'   # uni blue
    color2 = '#D20537'  # uni red
    uni_gray= '#BEC3C8'
    ax1.set_xlabel(x_label)
    ax1.set_ylabel('Average interligand distance (Å)')
    ax1.plot(plot_x,plot_y_plus_std, color='white')
    ax1.plot(plot_x,plot_y_minus_std, color='white')
    ax1.plot(plot_x,plot_y,color=color, linewidth=3)
    plt.fill_between(plot_x,plot_y_plus_std,plot_y_minus_std, color=color, alpha=0.2)

    if times_inclusion == True:
        ax1.set_xticks(x_ticks)

    if plot_limits == False:
        x_limit_lower = plot_x[0]
        x_upper_limit = plot_x[-1]

        if times_inclusion == True:
            x_limit_upper = max_simulation_time     # make sure last timepoint is last timepoint of sim

        y_limit_lower = 0
        y_limit_upper = 120

        plt.xlim(x_limit_lower, x_limit_upper)         # setting plot limits for x axis
        plt.ylim(y_limit_lower, y_limit_upper)

    else:
        pass

    # Saving figure of plot
    plt.savefig(working_directory + '_' + identifier + '_resultplot.png', dpi=300)
    print('Saved figure. Done.')
