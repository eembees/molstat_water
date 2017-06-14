import openbabel as ob
import pybel as pb
import numpy as np

global mol
global FF


def get_energy(filename):
    """ Returns the MMFF94 energy of the molecule in kcal/mol.
    """

    global mol

    mol = pb.readfile('xyz', filename).next()

    FF = ob.OBForceField.FindForceField("MMFF94")
    FF.Setup(mol.OBMol)

    return FF.Energy()

"""
def generate_chain(n):
""" #Generate a n-length carbon chain.
"""

    global mol

    carbon = "C"*n
    mol = pb.readstring("smi", carbon)
    mol.addh()
    mol.make3D()

"""
def set_dihedral(angles):
 #Set the dihedral angles of the carbon chain


    global mol
    global forcefield

    constraints = ob.OBFFConstraints()
    # constraints.AddDistanceConstraint(1, 10, 3.4)       # Angstroms
    # constraints.AddAngleConstraint(1, 2, 3, 120.0)      # Degrees
    # constraints.AddTorsionConstraint(1, 2, 3, 4, 180.0) # Degrees


    # Find all carbons
    smarts = pb.Smarts("C")
    smarts = smarts.findall(mol)

    for i in xrange(len(angles)):

        angle = angles[i]

        # Get the next 4 carbons
        Cs = smarts[i:i+4]

        ai, = Cs[0]
        bi, = Cs[1]
        ci, = Cs[2]
        di, = Cs[3]

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


def find_local_min(steps):
    """
    """

    global mol
    global FF

    FF = ob.OBForceField.FindForceField("MMFF94")
    # tested both 5 and 25 steps. As far as I remember, 5 steps were enough
    FF.SteepestDescent(steps)
    FF.GetCoordinates(mol.OBMol)


def save_molecule(filename):
    """ Save the molecule as <filename> in XYZ format.
    """

    global mol

    mol.write('xyz', filename, overwrite="True")


if __name__ == '__main__':

    generate_chain(6)

    print get_energy()

    # save_molecule('carbon_man1.xyz')

    set_dihedral([180.0, -180.0, 180.0])

    print get_energy()

    # save_molecule('carbon_man2.xyz')

    find_local_min()

    print get_energy()

    # save_molecule('carbon_man3.xyz')
