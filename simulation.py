## simulation.py Version 1.0
## Written by Magnus Berg Sletfjerding (eembees) and
###########################################################
"""Importing modules"""
import openbabel as ob
import pybel as pb
import numpy as np
import filecmp
import copy
import os.path
from math import pi, sin, cos

import extract as ex
import rotation as rot
import movement as mov
import molec as mc

"""Running Script for simulation"""
# # Defining universal variables
# initial_energy = mc.get_energy()
n_steps = 50
energies_before = []
energies_after = []
mol_name, elements, coordinates = ex.readfile('w6.xyz')
molecules = ex.divide(elements, coordinates)
d = 1
writing = 1 # # CHANGE THIS TO 1 TO WRITE FILES



"""
First step: Creating series of randomly generated files
"""

for i in range(n_steps):
    # # Generating random value, deciding on rotation or movement
    rot_not = np.random.choice([0,1])

    if rot_not == 1:
        """
        Rotates molecules by the same random angle.
        """
        n_molecules_rotating = np.random.randint(1,6)
        angle = 2 * pi * np.random.rand()           # Consider randomizing angle completely
        for j in range(n_molecules_rotating):
            # # Defining which molecule rotates
            rot_mol_num = np.random.choice(6)
            rot_mol_cor = molecules[rot_mol_num]

            # # Defining which atom rotates, and which are axes
            rot_atom_num = np.random.choice(3)
            rot_atom = rot_mol_cor[rot_atom_num]
            axis_1 = rot_mol_cor[(rot_atom_num + 1) % 3]
            axis_2 = rot_mol_cor[(rot_atom_num + 2) % 3]

            # # Executing rotation3d
            rot_atom_done = rot.rotation3d(axis_1, axis_2, rot_atom, angle)
            # # Restoring and rewriting coordinate of rotated atom
            molecules[rot_mol_num][rot_atom_num] = rot_atom_done

            write_now = 1

    elif rot_not == 0:
        """
        Moves Molecules
        """
        n_molecules_moving = np.random.randint(1,6)
        for j in range(n_molecules_moving):
            # # Defining rotating molecule
            mov_mol_num = np.random.choice(6)
            mov_mol_old = copy.copy(molecules[mov_mol_num])

            # # Executing movement
            mov_mol_new = mov.randommove(molecules[mov_mol_num], d)

            # print mov_mol_new
            # # Limiting movement
            if all(x < 0.0 for x in (mov_mol_new[0][0],
                                     mov_mol_new[1][0],
                                     mov_mol_new[2][0])) and \
               all(x > -7.0 for x in (mov_mol_new[0][0],
                                     mov_mol_new[1][0],
                                     mov_mol_new[2][0])) and \
               all(x < 3.0 for x in (mov_mol_new[0][1],
                                     mov_mol_new[1][1],
                                     mov_mol_new[2][1])) and \
               all(x > -3.0 for x in (mov_mol_new[0][1],
                                      mov_mol_new[1][1],
                                      mov_mol_new[2][1])) and \
               all(x < 3.5 for x in (mov_mol_new[0][0],
                                     mov_mol_new[1][0],
                                     mov_mol_new[2][0])) and \
               all(x > -3.0 for x in (mov_mol_new[0][0],
                                      mov_mol_new[1][0],
                                      mov_mol_new[2][0])):

                # Accept change:
                molecules[mov_mol_num] = mov_mol_new
                write_now = 1
            else:
                write_now = 0
            # molecules[mov_mol_num]


    # # Writing to new file
    if writing == 1 and write_now == 1:
        newfilename = "w6_%s" % i + ".xyz"
        coordinates = ex.unite(molecules)
        ex.writefile(newfilename, mol_name, elements, coordinates) # TODO: find how to configure inputs correctly, with n_steps as part of filename


"""
PART II: optimizing energy of the xyz files
"""
    # # Getting energy and optimizing
for i in range(n_steps):
    filename = "w6_%s" % i + ".xyz"
    if os.path.isfile('./%s'%filename):
        energies_before.append(mc.get_energy(filename))

        mc.find_local_min()

        energies_after.append(mc.get_energy(filename))
        mc.save_molecule(filename)

print 'energies after are \n', np.asarray(energies_after), min(np.asarray(energies_after))
