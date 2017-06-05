"""This code is inspired by molecule.py by Jan Jensen, as well as code found at
'https://stackoverflow.com/questions/17763655/rotation-of-a-point-in-3d-about-an-arbitrary-axis-using-python'
Importing modules"""
import openbabel as ob
import pybel as pb
import numpy as np
from math import pi, sin, cos

global mol
global forcefield

"""Defining Functions"""

def set_dihedral(angles):
    """ Modifying set_dihedral from Genetic Algoritm project, simulating alkanes
        Setting the dihedral angles of the molecule, using input angles list
    """

    global mol
    global forcefield

    constraints = ob.OBFFConstraints()
    # constraints.AddDistanceConstraint(1, 10, 3.4)       # Angstroms
    # constraints.AddAngleConstraint(1, 2, 3, 120.0)      # Degrees
    # constraints.AddTorsionConstraint(1, 2, 3, 4, 180.0) # Degrees


    # Find all Oxygen atoms
    smarts = pb.Smarts("O")
    smarts = smarts.findall(mol)

    for i in xrange(len(angles)):

        angle = angles[i]

        # Get the next 4 Oxygens
        Os = smarts[i:i+4]

        ai, = Os[0]
        bi, = Os[1]
        ci, = Cs[2]
        di, = Os[3]

        a = mol.OBMol.GetAtom(ai)
        b = mol.OBMol.GetAtom(bi)
        c = mol.OBMol.GetAtom(ci)
        d = mol.OBMol.GetAtom(di)

        mol.OBMol.SetTorsion(a, b, c, d, angle/(180.0)*np.pi) # Radians
        anglep = mol.OBMol.GetTorsion(a, b, c, d)

        # Define constraint
        constraints.AddTorsionConstraint(ai, bi, ci, di, anglep) # Degrees

    # Setup the force field with the constraints
    forcefield = ob.OBForceField.FindForceField("GAFF")
    forcefield.Setup(mol.OBMol)
    forcefield.SetConstraints(constraints)
    forcefield.EnableCutOff(True) # VDW cutoff
    forcefield.SetElectrostaticCutOff(0) # Remove electrostatics

    # Use forcefield to reoptimize bondlengths+angles
    find_local_min()


"""Rotation Functions"""


def R(theta, u):
    """Defines the 3D rotation matrix for rotaion about unit vector passing through origin.
    u is unit vector, with 3 elements
    theta is angle of Rotation"""
    return [[cos(theta) + u[0]**2 * (1-cos(theta)),
             u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta),
             u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
            [u[0] * u[1] * (1-cos(theta)) + u[2] * sin(theta),
             cos(theta) + u[1]**2 * (1-cos(theta)),
             u[1] * u[2] * (1 - cos(theta)) - u[0] * sin(theta)],
            [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
             u[1] * u[2] * (1-cos(theta)) + u[0] * sin(theta),
             cos(theta) + u[2]**2 * (1-cos(theta))]]

def Rotate(pointToRotate, point1, point2, theta):

    """INPUTS:  point1 and point2: coordinates of the two other atoms in the molecule
       OUTPUTS: new position of point"""

    u= []
    squaredSum = 0
    for i,f in zip(point1, point2):
        u.append(f-i)
        squaredSum += (f-i) **2

    u = [i/squaredSum for i in u]

    r = R(theta, u)
    rotated = []

    for i in range(3):
        rotated.append(round(sum([r[j][i] * pointToRotate[j] for j in range(3)])))

    return rotated
