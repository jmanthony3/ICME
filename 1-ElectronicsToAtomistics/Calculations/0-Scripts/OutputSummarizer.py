#!/usr/bin/env python

import numpy as np
import os
import re
from time import time
from subprocess import Popen, PIPE


#
### EDIT THIS PARAMETER FOR CODE IN USE
# 'qe' - Quantum Espresso
# 'vasp' - VASP
dft = 'qe'

# Number of processors
num_proc = 16

##### MATERIAL SETTINGS - USER EDIT REQUIRED #####

# vvv These parameters are for QE only, VASP will use the available POTCAR file vvv
# Directory containing your pseudopotentials
pp_dir = '.'
# Element name
el = 'Fe'
# Potential file name
potential = 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF'
# Element weight
el_weight = 55.845
# ^^^


##### DFT PARAMETERS - USER EDIT IF DESIRED #####
# Energy cutoff value (eV) (qe only)
energy_cutoff = 64*13.6057
# The below energy cutoff setting works for most materials.
# Check your pseudopotential file to see the suggested value for your material
energy_cutoff_rho = energy_cutoff * 4.0
# K-points specification
kpoints = 8
# Smearing parameter (QE only)
smear = 0.06
# Relax the ions in the z direction? (vasp only right now)
relax = False


##### STACKING FAULT PARAMETERS - USER EDIT IF NEEDED #####
# Number or location of simulated points. It MUST start with 0.
#fault_points = [0.0, 1./4., 1./2.]
fault_points = np.linspace(0,1,16)
# Size of vacuum in Ang.
vacuum = 20.0
# Number of stacking layers - CAREFUL, DFT scales very poorly with more atoms
stacking_layers = 10


##### Unit conversions #####
ry_to_ev = 13.6056849587
au_to_ang = 0.52917721092

##### Below this line should not require user edits #####

# Function that creates the geometry and writes the input file
def create_stackingfault( struct, lp, layers, slip, system=None ):
	# struct : reference structure
	# lp : lattice parameter\
	# layers : number of layers of atoms to generate
	# slip : normalized displacement in slip direction
	# system : the type of stacking fault to create

	if struct.lower() == 'fcc':
		basis = np.array( [[ np.sqrt(2.)/4., np.sqrt(6.)/4., 0.0 ],
				   [ -np.sqrt(2.)/4., np.sqrt(6.)/4., 0.0 ],
				   [ 0.0, np.sqrt(6.)/6., np.sqrt(3.)/3.]] )
		base_atoms = np.array( [[ 0.0, 0.0, 0.0]] )
		if not system or system == 'partial' or system == 'shockley':
			# Assume the partial dislocation is desired (most common)
			slip_vector = np.array([0.0, np.sqrt(6.)/2., 0.0]) # Partial dislocation line
		elif system == 'full' or system == 'burgers':
			slip_vector = basis[0,:] # burgers vector direction

	elif struct.lower() == 'bcc':
		basis = np.array( [[ (np.sqrt(2.)/4. + 0.5), -(0.5 - np.sqrt(2.)/4.), 0.0 ],
				   [ -(0.5 - np.sqrt(2.)/4.), (np.sqrt(2.)/4. + 0.5), 0.0 ],
				   [ np.sqrt(2.)/4., np.sqrt(2.)/4., np.sqrt(2.)/2.]] )
		base_atoms = np.array( [[ 0.0, 0.0, 0.0]] )
		if not system or system == 'partial':
			# Assume the partial dislocation is desired (most common)
			slip_vector = np.array([np.sqrt(2.)/2.,np.sqrt(2.)/2., 0.0]) # Partial dislocation line?
		elif system == 'longpartial':
			# or?
			slip_vector = np.array([0.5,-0.5,0.0])
		elif system == 'full' or system == 'burgers':
			slip_vector = basis[0,:] # burgers vector directionn

	# Create atoms
	atoms = []
	threshold = np.cross(basis[0,0:2],basis[1,0:2])
	if threshold < 0:
		threshold = np.cross(basis[1,0:2],basis[0,0:2])
		d1 = lambda x : -basis[1,1]*x[0] + basis[1,0]*x[1]
		d2 = lambda x : -basis[0,1]*x[0] + basis[0,0]*x[1]
		b1 = basis[1,:]
		b2 = basis[0,:]
	else:
		d1 = lambda x : -basis[0,1]*x[0] + basis[0,0]*x[1]
		d2 = lambda x : -basis[1,1]*x[0] + basis[1,0]*x[1]
		b1 = basis[0,:]
		b2 = basis[1,:]

	for i in range(layers):
		if i == 0:
			old_atom = np.zeros(3)
			atoms.append(old_atom)
			continue

		new_atom = old_atom + basis[2,:]

		# Check if atom is within cell boundaries
		in_cellx, in_celly = False, False
		while not in_cellx or not in_celly:
			if d1(new_atom) < 0:
				new_atom += b2
			elif d1(new_atom) > threshold:
				new_atom -= b2
			else:
				in_cellx = True

			if d2(new_atom) > 0:
				new_atom += b1
			elif d2(new_atom) < -threshold:
				new_atom -= b1
			else:
				in_celly = True

		old_atom = new_atom.copy()

		if i >= int(layers/2):
			new_atom += slip*slip_vector

			# Check if atom is within cell boundaries
			in_cellx, in_celly = False, False
			while not in_cellx or not in_celly:
				if d1(new_atom) < 0:
					new_atom += b2
				elif d1(new_atom) > threshold:
					new_atom -= b2
				else:
					in_cellx = True

				if d2(new_atom) > 0:
					new_atom += b1
				elif d2(new_atom) < -threshold:
					new_atom -= b1
				else:
					in_celly = True

		atoms.append(new_atom)

	if dft == 'vasp':
		write_vasp_inputs(lp, basis, layers, atoms)
	elif dft == 'qe':
		write_qe_inputs(lp, basis, layers, atoms)

	# Returns the stacking fault area
	return [np.linalg.norm(np.cross(basis[0,:],basis[1,:])), np.linalg.norm(slip_vector), len(atoms)]
# end create_stackingfault

def write_qe_inputs(lp, basis, layers, atoms):

	# Write the input file
	with open('gsfe.in','w') as f:
		# Control section
		f.write(' &control' + os.linesep)
		f.write("\tprefix=''" + os.linesep)
		f.write("\toutdir='temp'" + os.linesep)
		f.write("\tpseudo_dir = '{}',".format(pp_dir) + os.linesep)
		f.write(' /' + os.linesep)
		# System section
		f.write(' &system' + os.linesep)
		f.write('\tibrav= 0, nat= {}, ntyp= 1,'.format(layers) + os.linesep)
		f.write('\tcelldm(1) ={0}, '.format(lp/au_to_ang) + os.linesep)  # Convert lattice parameter to a.u.
		f.write('\tecutwfc ={},ecutrho ={}, '.format(energy_cutoff/ry_to_ev,\
		 	energy_cutoff_rho/ry_to_ev) + os.linesep)  # Convert energy to Rydberg
		f.write("\toccupations='smearing', smearing='mp', degauss={} ".format(smear)\
		 	+ os.linesep)
		f.write(' / ' + os.linesep)
		# Electrons section
		f.write(' &electrons ' + os.linesep)
		f.write("mixing_mode ='local-TF', electron_maxstep = 110," + os.linesep)
		f.write(' / ' + os.linesep)
		# Atomic species - pseudopotential specification
		f.write('ATOMIC_SPECIES ' + os.linesep)
		f.write(' {}  {} {} '.format(el,el_weight,potential) + os.linesep)
		# Atomic positions
		f.write('ATOMIC_POSITIONS alat ' + os.linesep)
		for a in atoms:
			f.write(' {1}\t{0[0]}\t{0[1]}\t{0[2]} '.format(a, el) + os.linesep)
		# K-points
		f.write('K_POINTS automatic ' + os.linesep)
		f.write(' {0} {0} 1 0 0 0 '.format(kpoints) + os.linesep)
		# Basis vectors
		f.write('CELL_PARAMETERS alat ' + os.linesep)
		f.write('{0[0]}\t{0[1]}\t{0[2]} '.format(basis[0,:]) + os.linesep)
		f.write('{0[0]}\t{0[1]}\t{0[2]} '.format(basis[1,:]) + os.linesep)
		f.write('0.0\t0.0\t{} '.format(basis[2,2]*10. + vacuum/lp) + os.linesep) # Adds vacuum.
# end write_qe_inputs

def write_vasp_inputs(lp, basis, layers, atoms):

	# Write the POSCAR file
	with open('POSCAR','w') as f:
		f.write(el + os.linesep)
		f.write('%f '%lp + os.linesep)
		f.write('{0[0]}\t{0[1]}\t0.0 '.format(basis[0,:]) + os.linesep)
		f.write('{0[0]}\t{0[1]}\t0.0 '.format(basis[1,:]) + os.linesep)
		f.write('0.0\t0.0\t{0} '.format(basis[2,2]*layers + 20.0/lp) + os.linesep)
		f.write('%i '%layers + os.linesep)
		f.write('Cartesian ' + os.linesep)
		for a in atoms:
			f.write('{0[0]}\t{0[1]}\t{0[2]}\tF\tF\tT '.format(a) + os.linesep)

	# Write the KPOINTS file
	with open('KPOINTS', 'w') as f:
		f.write(el + os.linesep)
		f.write('0 ' + os.linesep)
		f.write('Monkhorst-Pack ' + os.linesep)
		f.write('{0} {0} 1 '.format(kpoints) + os.linesep)
		f.write('0 0 0 ' + os.linesep)
# end write_vasp_inputs

def gsfe( struct, lp, slip_system=None ):

	# Creates the summary file
	with open('GSFE_SUMMARY','w') as f:
		pass
		#f.write("Generalized Stacking Fault Energy for {} {}\n".format(struct,el))
		#f.write("=============================================\n")

	# Initialize the lists for the energy values
	energy = []
	filenum = 1
	for d in fault_points:
		# Create the stacking fault structure
		[area, fault_length, natoms] = create_stackingfault( struct, lp, stacking_layers, d, slip_system )

		if dft == 'qe':
			E, walltime = run_qe(filenum)
			E = E*ry_to_ev
		

		energy.append(E)

		with open('GSFE_SUMMARY','a') as f:
			# Write the normalized displacement and the relaxed and static energies in mJ/m**2
			f.write('{}\t{}\t{} '.format(d*fault_length,
					(energy[-1] - energy[0]) * 1.60217733e-19 * 1e23 / area,
					walltime) + os.linesep)
		filenum = filenum + 1
# end gsfe

def run_qe(filenum):

		# Run quantum espresso
		# os.system('mpirun -np {} pw.x < gsfe.in > gsfe.out'.format(num_proc))
		# print("Opening 'gsfe_{}.out'".format(filenum))

		# Get energy
		p1 = Popen(['grep','! *[ ] total energy',"gsfe_{}.out".format(filenum)], stdout=PIPE)
		p2 = Popen(['awk','{print $5}'], stdin=p1.stdout, stdout=PIPE)
		try:
			foo = float(Popen(['awk','{print $5}'], stdin=Popen(['grep','! *[ ] total energy',"gsfe_{}.out".format(filenum)], stdout=PIPE).stdout, stdout=PIPE).communicate()[0])
		except ValueError:
			foo = float(Popen("sed -n 's%\!*[[:space:]]*total energy[[:space:]]* = [[:space:]]*%%p' gsfe_{}.out | sed 's% Ry$%%' | tail -1".format(filenum), shell=True, stdout=PIPE).stdout.read())

		# Get walltime
		try:
			p3 = Popen(['grep','PWSCF',"gsfe_{}.out".format(filenum)],stdout=PIPE)
			p4 = Popen(['tail','-n1'],stdin=p3.stdout,stdout=PIPE)

			time_str = p4.communicate()[0]
			time_arr = re.findall('(\d+)m[ ]*(\d+).(\d+)s', time_str)[-1]
			time_arr = [int(i) for i in time_arr]

			if len(time_arr) == 3:
				time = 60*time_arr[0] + time_arr[1] + 0.01*time_arr[2]
			else:
				time = time_arr[0] + 0.01*time_arr[1]
		except:
			time = 0.0

		# print([foo*ry_to_ev, time])
		# print("Closing gsfe_{}.out'".format(filenum))
		return [foo, time]
# end run_qe

def run_vasp():

		# Run vasp - Static run first.
		tries = 0
		has_run = False
		while tries < 2 and not has_run:
			try:
				if relax:
					# Run vasp - first a relaxation, then a static run
					os.system('cp relax.INCAR INCAR')
					msg = os.system('mpirun -np {} ./vasp'.format(num_proc))
					os.system('cp CONTCAR POSCAR')
				os.system('rm WAVECAR CHGCAR')
				os.system('cp gsfe.INCAR INCAR')
				msg = os.system('mpirun -np {} ./vasp'.format(num_proc))
				p1 = Popen(['tail','-n1','OSZICAR'], stdout=PIPE)
				p2 = Popen(['awk','{print $5}'], stdin=p1.stdout, stdout=PIPE)
				E = float(p2.communicate()[0])
				has_run = True
			except Exception:
				print("VASP run failed, trying again after deleting output files.")
				try:
					os.remove('IBZKPT CHG CONTCAR DOSCAR EIGENVAL OSZICAR OUTCAR'\
				+ ' PCDAT XDATCAR EIGENVAL vasprun.xml cellvol PROCAR CHGCAR WAVECAR')
				except Exception:
					pass
				tries += 1
		# Get energy
		p1 = Popen(['tail','-n1','OSZICAR'], stdout=PIPE)
		p2 = Popen(['awk','{print $5}'], stdin=p1.stdout, stdout=PIPE)

		return float(p2.communicate()[0])


		#if msg != 0:
		#	print("VASP exited with an error. It may only be able to run on the CAVS cluster machines.")
# end run_vasp


if __name__ == "__main__":
	import sys
	argc = len(sys.argv)

	# Get inputs from the command line
	struct = sys.argv[1] if argc > 1 else None
	latp = float(sys.argv[2]) if argc > 2 else None
	slip = sys.argv[3] if argc > 3 else None
	extent = float(sys.argv[4]) if argc > 4 else None

	# If the inputs are not specified, prompt for them.
	#if not element:
	#	print("Enter the element name:")
	#	element = raw_input("> ")
	if not struct:
		i_struct = None
		while not i_struct:
			print("== Enter the desired structure ==")
			print("= 1) FCC                        =")
			print("= 2) BCC                        =")
			print("=================================")
			i_struct = int(raw_input("> "))
			if i_struct == 1:
				struct = 'fcc'
				print("Structure type: FCC")
			elif i_struct == 2:
				struct = 'bcc'
				print("Structure type: BCC")
			else:
				print(" Only FCC and BCC systems are implemented currently ")
				i_struct = None
	if struct.lower() == 'fcc':
		slip_choices = ['(111)[1-10] (full)','(111)[11-2] (partial)']
	elif struct.lower() == 'bcc':
		slip_choices = ['(110)[-111] (full)','(110)[001] (partial)','(110)[-110] (longpartial)']
	else:
		print("Structures other than FCC and BCC are not currently supported")
		sys.exit(1)
	if not latp:
		print("Enter the equilibrium lattice parameter for your element:")
		latp = float(raw_input("> "))
	if not slip:
		i_slip = None
		while not i_slip:
			print("====== Enter the desired slip system ======")
			for i,sc in enumerate(slip_choices):
				print("= {0}) {1:<36} =".format(i+1,sc))
			print("===========================================")
			i_slip = int(raw_input("> "))
			if i_slip == 1:
				slip = 'full'
				print("Using burgers vector direction")
			elif i_slip == 2:
				slip = 'partial'
				print("Using partial dislocation direction")
			elif i_slip == 3 and struct == 'bcc':
				slip = 'longpartial'
				print("Using longer partial dislocation direction")
			else:
				print(os.linesep + "Choice not recognized" + os.linesep)
				i_slip = None
	if not extent:
		extent = 1.0

	gsfe(struct, latp, slip)
