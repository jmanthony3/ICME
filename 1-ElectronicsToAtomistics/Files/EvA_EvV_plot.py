# The purpose of this code is to create plots for any set of EvsA and EvsV data
# It should be ran from a directory that has EvsA and EvsV data in it
# This code can both print the figures to the screen and write the files to .pdf 
# depedning on user settings

import numpy as np
import matplotlib.pyplot as plt
import os
# The values are stored in 2 Columns as:
# Lattice Parameter(angstrom)       energy(eV)
path = os.getcwd()



#######################################################
###                  User Settings                  ###
#######################################################

# Print to screen/pdf or not
write_to_screen = True
print_to_pdf    = True 

energy_offset = 2858.8298734 # Set by running simulation with very large lattice parameter

#Color options for plots
# Best source for colors = https://xkcd.com/color/rgb/
color1 = '#b66325'  
color2 = '#8b2e16'  
color3 = '#00035b'  

mk_size = 8  # Markersize for plots, larger is bigger
ax_label_size   = 36 # Font size for axis labels
title_font_size = 42 # Font size for table title
marker_type = 'o'  # Type of marker used for plotting google matplotlib for more options

EvA_out_name = "Name_of_EvA.pdf"
EvV_out_name = "Name_of_EvV.pdf"
Combined_out_name = "Name_of_Combined.pdf"

#Path to files.  the / is for linux. If you are on windows change them to \
EsvA = open(path+'/EvsA') #should point to EvsA.txt file
EsvV = open(path+'/EvsV') #should point to EvsV.txt file


#######################################################
###            End of  User Settings                ###
#######################################################


#Initializing arrays
A  = [] # array of Lattice Parameter in Angrstrom
A3 = [] # array of Volume in Angrstrom^3
E  = [] # array of Energy in eV



line = True
while line:  # looping over all of the lines in the file
    line = EsvA.readline()
    values = line.split()
    if not values:
        break
    A.append(float(values[0]))
    E.append(float(values[1]))
EsvA.close()


line = True
while line:   # looping over all of the lines in the file
    line = EsvV.readline()
    values = line.split()
    if not values:
        break
    A3.append(float(values[0]))
EsvV.close()


latt_param = np.array(A)        # creates array of values for lattice parameter
energy     = np.array(E)        # creates array of values for lattice parameter
energy = energy + energy_offset # shifting the energy values
volume = np.array(A3)


#creating plots

EvA, ax = plt.subplots(figsize = (12,10))
ax.plot(latt_param, energy, marker_type, markersize = mk_size, color = color1)
ax.axhline(linewidth=1, color='k')
ax.set_xlabel('Lattice Parameter ($\AA$)', fontsize = ax_label_size)
ax.set_ylabel('Energy ($eV$)' , fontsize = ax_label_size)
ax.set_title('Energy vs Lattice Parameter \n Curve for FCC copper', fontsize = title_font_size)
ax.grid(True, linestyle='-.')
ax.tick_params(labelcolor='k', labelsize='large', width=3)


EvV, ax = plt.subplots(figsize = (12,10))
ax.plot(volume, energy, marker_type, markersize = mk_size, color = color1)
ax.axhline(linewidth=1, color='k')
ax.set_xlabel('Volume ($\AA^{3}$)', fontsize = ax_label_size)
ax.set_ylabel('Energy ($eV$)' , fontsize = ax_label_size)
ax.set_title('Energy vs Volume \n Curve for FCC Copper', fontsize = title_font_size)
ax.grid(True, linestyle='-.')
ax.tick_params(labelcolor='k', labelsize='large', width=3)



combined , axarr = plt.subplots(1,2, figsize = (20,8))
combined.suptitle('Energy vs Lattice Parameter and Energy vs Cell Volume \n', fontsize = ax_label_size+4)
axarr[0].plot(latt_param, energy, marker_type, markersize = mk_size, color = color2)
axarr[0].set_xlabel('Lattice Parameter ($\AA$)', fontsize = ax_label_size)
axarr[0].set_ylabel('Energy (eV)' , fontsize = ax_label_size)
axarr[0].grid(True, linestyle='-.')
axarr[0].tick_params(labelcolor='k', labelsize=15, width=3)

axarr[1].plot(volume, energy, marker_type,markersize = mk_size, color = color2)
axarr[1].set_xlabel('Cell Volume ($\AA^{3}$)', fontsize = ax_label_size)
axarr[1].grid(True, linestyle='-.')
axarr[1].tick_params(labelcolor='k', labelsize=15, width=3)


# Printing plots to a .pdf file
if print_to_pdf == True:
    EvV.savefig(EvV_out_name, bbox_inches = 'tight')
    EvA.savefig(EvA_out_name, bbox_inches = 'tight')
    combined.savefig(Combined_out_name, bbox_inches = 'tight')

# Print to screen
if write_to_screen == True:
    plt.show()

