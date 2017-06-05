"""FUNCTIONS:
readfile -- Extracting values from .xyz file into arrays
writefile-- Writing the the file from coordinates arrays
divide   -- Categorizing array into subarrays
"""
# # Import modules
import numpy as np

def readfile(filename):
    """Reads the file to extract xyz data
    INPUTS: filename
    OUTPUTS: *name of molecule inxyz file (line 1)
             *names of elements, first column in xyz file
             *coordinates as ndarray
    """
    read = open(filename)
    lines = read.readlines()
    n_atoms = int(lines[0])
    mol_name = lines[1].strip()

    coordinates = np.zeros([n_atoms, 3])    # # Makes empty numpy array for all xyz coordinates
    elements = []                           # # Array listing all element symbols

    for i in range(len(coordinates)):
        # # extract the right line from file
        current_line = lines[i+2]
        current_line = current_line.split()

        # # write element names to element list
        elements.append(current_line[0])

        # # Write coordinates into ndarray
        current_coordinate = coordinates[i]
        for j in range(3):
            current_coordinate[j] = float(current_line[j+1])

        coordinates[i] = current_coordinate

    return mol_name, elements, coordinates

def writefile(filename, mol_name, elements, coordinates):
    fil = open(filename, 'w')
    fil.write("%d\n%s\n" % (coordinates.size / 3, mol_name))

    for i in range(len(elements)):
        current_coordinate = coordinates[i]
        fil.write("%s         %8.5f       %8.5f       %8.5f\n" % (elements[i], current_coordinate[0], current_coordinate[1], current_coordinate[2]))
    pass

def divide(elements, coordinates):
    """
    Divides coordinate list evenly among molecules
    """

    # Define uni. values
    n_atoms = int(len(coordinates))
    n_molecules = elements.count('O')
    if n_atoms % n_molecules != 0:
        print "ERROR: molecules and atoms do not match"
        print "Refer to line %s in extract.py to correct" %inspect.currentframe().f_back.f_lineno
        pass
    elements_per_molecule = n_atoms / n_molecules

    # Split array to smaller arrays
    molecules = []
    for i in range(1, n_molecules + 1):
        molecule = coordinates[(i-1)*elements_per_molecule:i*elements_per_molecule]
        molecules.append(molecule)

    return np.asarray(molecules)

def unite(molecules):
    """ 
    Unites coordinate list given molecule array as output from divide function
    Verified that united out array is identical to divide input array
    """
    n_molecules = len(molecules)
    coordinates = []
    for i in range(n_molecules):
        for j in range(len(molecules[0])):
            coordinates.append(molecules[i][j])
    return np.asarray(coordinates)
