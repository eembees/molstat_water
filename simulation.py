## simulation.py Version 1.0
## Written by Magnus Berg Sletfjerding (eembees) and
###########################################################
"""Importing modules"""
import openbabel as ob
import pybel as pb
import numpy as np
import filecmp
from math import pi, sin, cos

import extract as ex
import rotation as rot
import movement as mov
import molecule as mc


"""Running Script for simulation"""
# # Defining universal variables
initial_energy = mc.get_energy('w6.xyz')
n_steps = 100
energies = []
mol_name, elements, coordinates = readfile('w6.xyz')
molecules = ex.divide(elements, coordinates)


"""
First step: Creating series of randomly generated files
"""

for i in n_steps:
    # # Generating random value, deciding on rotation or movement
    rot_not = np.random.getrandbits(1)

    if rot_not == 1:
        # Rotates molecules
        n_molecules_rotating = np.random.randint(1,6)
        angle = pi * np.ra
        for j in range(n_molecules_rotating):
            rot_mol = np.random.choice(6)
            rot_mol_cor = molecules[rot_mol]
            rot_mol_cor



    elif rot_not == 0:
        # Moves Molecules
        n_molecules_moving = np.random.randint(1,6)
        for j in range(n_molecules_moving):
            mov_mol = np.random.choice(6)


    # # Writing to new file
    newfilename = "w6_%s" % i + ".xyz"
    newcoordinates = ex.unite(molecules)
    ex.writefile(newfilename, mol_name, elements, coordinates) # TODO: find how to configure inputs correctly, with n_steps as part of filename

    # # Getting energy and optimizing



    energies.append
