[//]: # 'Report for molstat'

# Acknowledgements
My thanks go to Jan Jensen, for suggestions and guidance throughout the project, Freja Kellerup [?] for the necessary introduction to the different python modules we used, and of course to Michelle Berling, without whom this project would never be completed in time.  

# Introduction to the water project
The overall water project aims to simulate the behavior of water molecules with each other, with the long term approach to involve different solutes, observing the interaction between molecules on the atomic level.
## Importance of the project
The project\'s main importance lies in the potential to develop a new methodology for simulations of water molecules. Finding a most accurate way to simulate the water molecule\'s movements is imperative to finding a holistic model, and thus developing, in the future, a more accurate methodology for simulating molecules in solution accurately. On the long term, this project can have impacts within every field that relies on simulations of molecules, including, but not limited to drug design and industrial catalyses.
# Theory

The success of the project, in essence, depends on the integrity of the energy calculations - without which the project would ultimately be moot. In line with the instructions, we used the MMFF94 force field to analyze the energies of the states of the molecules.

## MMFF94
The MMFF94 force field, as described by Halgren [CITE], takes seven parameters for the energy of a system, not of which all are relevant to our case, but nevertheless, the parameters include:
+ Bond stretching
+ Angle Bending
+ Stretch-bend interactions
+ Out-of-Plane Bending at Tricoordinate centers
+ Torsion interactions
+ Van Der Waals interactions
+ Electrostatic interactions
and are all part of the equation following:
[INSERT EQUATION 1]


The inclusion of all of these parameters in the report would be unneccesary, so the focus will be on just a few.
### Bond stretching
[INSERT EQUATION 2]
Here, kb\_ij references the force constant, while delta\_r\_ij references the deviance from the reference bond length in question. In short terms, when the bond length is equal to the reference, EB\_ij = 0.

### Van der Waals interactions

The MMFF94 field uses the "buffered-14-7" form, as described in 1992 [ Halgren J. Am. Chem. Soc.]. The advantage of the MMFF94 field is the fact that it takes VDW interactions into account when atoms are more than 3 bonds apart, or as in our case, parts of different domains. This feature allows for analyzing multi-domain interactions, and in our case, accurately observing the energy of different states.

### Electrostatic interactions
[EQUATION]
where d is  0.05 Å and D is the dielectric constant. In our case, D = 1.0, as the interactions between water molecules are dependent only on each other, and the space between the atoms are vacuum. The q's represent partial atomic charges, dependent on [EQUATION 14] seen in as [MORE IN PART V] etc..

SOURCES FOR THIS SECTION:

"Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94" Halgren et al 96

"Merck molecular force field. II. Basis, form, scope, parameterization, and performance of MMFF94" Halgren et al 96

[ Halgren J. Am. Chem. Soc.].

## Rotation
In order to best facilitate the rotation script, only the hydrogens were chosen to be rotated around the axis of the rest of the molecule.
In order to do this, we researched 3-dimensional rotation, and found a guide by Paul Bourke [cite], outlining the process, as well as how to find the rotation matrix, given an axis direction vector.
The rotation matrix we defined was found to be:
[WRITE MATRIX HERE]
where [x, y, z] represents the direction vector of the axis, defined as the normalized difference of the two axis points. When we multiply any point by this rotation matrix, the point would rotate around the axis by the angle defined in the matrix.
The rotation, of course, requires the transformation of the molecule to the origin before rotation can occur - as the rotation takes place, around a normalized direction vector, that is defined from the origin. Hence, a reference axis position vector needed to be first subtracted from the rotation point vector, and then added to the rotated point vector to get the final position vector of the rotated point.

SOURCES:
http://paulbourke.net/geometry/rotate/

## Movement
The theory behind the movement script is quite simple. Rather than construct a set of vectors with which the molecules could move, the movement was made completely random, with the option to define the distance of the movement. Hence, the movement of any point is executed by adding the position vector of the point with the movement vector, either randomly, or quasi-randomly generated.

# Algorithm

The overall algorithm is divided into five main parts:  
The modules we import into our simulation:
+ extract.py, the script that extracts and writes the files defining the states
+ rotation.py, the script containing the function that rotates molecules
+ movement.py, the script containing the movement function
+ molec.py, a modified version of *molecule.py* (provided beforehand), which includes the the minimization function, a function to write the file after minimization, as well as the MMFF94 energy calculation function, using the openbabel module  

The script used to run the simulation:
+ simulation.py

## extract.py
The four functions included in the*extract.py* script are the following:
+ **readfile**
+ **writefile**
+ **divide**
+ **unite**

### readfile and writefile
Reads and writes files to and from \*.xyz files. When reading the file, readfile will convert the file into a 3d array of arrays, where len(coordinates) = n\_atoms, and len(coordinates[0]) = 3. We chose this approach because it made it very simple to reference the number of the atom in question to the coordinates when treating the arrays. Furthermore, we then had the possibility, while bug testing, to ensure that all three dimensions, separately, were identical. For example, at one point, during the debugging and testing of *movement.py*, we found that it did not change the z-coordinate of the molecules, by comparing the three different arrays with Boolean operations. The **readfile** function also extracts the element names into a separate array, so that it may be written into the file using **writefile,** but also to be used as reference in other functions, such as the **divide** function.

The **writefile** function takes a new filename, the molecule name, the element array, and the coordinate arrays as arguments, and writes the files into a new \*.xyz file. The function was written manually, as we extracted the values manually, and hence wished to be in full control of the writing as well.

### divide and unite
The necessity of the **divide** function (and subsequently, the **unite function**) arose during rotation. It became increasingly prevalent making lists of the molecules to be rotated was far too time-consuming. Hence, we decided to split the molecules into discrete series of arrays, introducing another dimension in our ndarray. The divide function takes as arguments the element array and coordinates array from the **readfile** output, and outputs an ndarray of dimensions (n\_molecules, elements\_per\_molecule, 3). The divide function decides the number of molecules by counting the number of Oxygen atoms in the elements array. This obviously lessens the generalization opportunities of the code, but it does however, leave the possibility to modify the code to any sets of triatomic molecules (with the obvious exception of trianions, ozone, and simulations of thiozone).

The **unite** function, very simply, reverses the divide function, returning an array identical (at least in dimensions) to the output of the **readfile** function.



The *extract.py* also includes an "if __name__ == '__main__':" test, using the filecmp module to test whether it can write files that are identical.

## rotation.py
This is, by far, the most complex module, due to the complex mathematical operations it involves. Even so, it includes only one function - **rotation3d** - and an "if __name__ == '__main__':" test.
The **rotation3d** function takes three points and an angle (in radians) as its arguments. Two of these points are axis points, one is the point to be rotated. From there, the rotation goes in three main steps:
### Translating and Initializing
The first step is to define the position vector ax, which serves as the representation of the point while axis1 is defined as the origin. Hence, the rotation operation can assume the axis intersects the origin, making our job a whole lot easier.  
[ADD IMAGE]
The first step also initializes the rotation point vector, which serves to represent the point once it has been rotated around the origin. Finally, the direction vector is defined and normalized by the following code (lines 32-34).
~~~~
N = map(sum, zip(axis2, axis1_neg)) # axis vector
sqsum = (sum(sqrt(N[i]**2) for i in range(3)))
direction = [N[i]/sqsum for i in range(3)]
~~~~

### Rotating
This part begins by simplifying the 3d rotation matrix factors used, in order to clean up the matrix, providing maximum readability and customizability. Another benefit is that the trigonometric equations need only be calculated once for each rotation, instead of three times per rotation.   
Once the factors are simplified, the matrix D is defined as described above, using this snippet of code:
~~~
# Matrix 'D'[3x3] 3D rotation matrix
d11 = tr*x**2 + co
d12 = tr*x*y - si*z
d13 = tr*x*z + si*y
d21 = tr*x*y + si*z
d22 = tr*y**2 + co
d23 = tr*y*z - si*x
d31 = tr*x*z - si*y
d32 = tr*y*z + si*x
d33 = tr*z**2 + co
~~~
Subsequently, the rotation point is defined by multiplying the ax point by the D matrix.
~~~
# # Define rot
rot[0] = d11 * ax[0] + d12 * ax[1] + d13 * ax[2]
rot[1] = d21 * ax[0] + d22 * ax[1] + d23 * ax[2]
rot[2] = d31 * ax[0] + d32 * ax[1] + d33 * ax[2]
~~~
While the matrix operation could easily have been done natively in python, we chose to manually program the matrix operations, as it did not seem to make the program slower.  


### Tranformation and output
Finally, the transformation of the rotated point back to the original area, simply done by adding the axis1 position vector, as:
~~~
ax = point - axis1
~~~
and thus, the new coordinates should be
~~~
newpoint = rot + axis1
~~~
The function finally outputs newpoint as an array of coordinates in the same format as point, axis1, and axis2.

### Limitations
One problematic point we ran into was that when rotating molecules, some hydrogen atoms would rotate such that their rotated coordinate would be localized not around its oxygen, but around its second hydrogen. We suspect that this is due to the fact that the function calculates a "negative" rotation vector when rotating the 2nd H molecule, rather than the 1st. We then had the option of modifying *rotation.py* or to modify our *simulation.py* and chose the latter, as it seemed the simpler option. This is the reason for the following snippet of code in *simulation.py*:  
~~~
if rot_atom_num == 1:
    axis_num_1 = 0
    axis_num_2 = 2
if rot_atom_num == 0:
    axis_num_1 = 2
    axis_num_2 = 1
if rot_atom_num == 2:
    axis_num_1 = 0
    axis_num_2 = 1
~~~
and this section seemed a more natural part to address this.

## movement.py
This module includes the **randommove** function.
**randommove** takes a set of points, no matter how large (it can be one atom, a molecule, multiple molecules, or the entire collection of molecules), and moves them randomly along the same random vector, of which the length can be varied. It takes the points, in an ndarray (as seen earlier, in the output of **readfile**), along with the requested distance in Å, as inputs. It does not include the normal if __name__ == __main__: test, as the function was shown to work on the first attempts, integrated in the *simulation.py* script.
The function begins by defining a direction vector, then finds a movemement vector by multiplying with the distance given, and adds the movement vector to the points. Next, it outputs the array newpoints, which is the same series of points, moved by the movement vector.

## molec.py
We modified the given *molecule.py* to our own *molec.py*, by introducing a few changes.
The first was to change the **get\_energy** function, which returns the energy of the state saved in an \*.xyz file.  
Second, we modified the **find\_local\_min** function to take the amount of steps as input, so it would be easier to change in the final simulation.
Although there were no more changes made, it is worth mentioning that the **find\_local\_min** function was always followed by the function **save_molecule**, ensuring that the state would always be saved. The reason we did not use openbabel/pybel in our transformations was simply that it seemed a better idea to do the rotation and movements ourselves - due to the fact that openbabel/pybel is more strictly defined, and we already had sufficient knowledge concerning vector operations and 3d rotations.

## simulation.py
The *simulation.py* script runs the simulation in a series of steps, as we need to read the files twice - once to manipulate the states, once to minimize, using the **find\_local\_min** function.

### Step 1: Generate random states
The script begins by a random decision whether to rotate or move molecules, using *numpy.random.choice* to find the variable rot\_not as either 1 or 0, which decides the outcome of the step.

#### Rotation
First, the rotation section decides on the number of molecules to rotate - this way, we can make a rotation transformation on up to all molecules simultaneously.
Second, the rotation defines the angle to rotate all molecules by, using another numpy random number.
Third, the script is setup to loop over the number of molecules to rotate, and randomly, for each iteration of the loop, chooses a molecule to rotate. It then chooses which H atom in the molecule to rotate (numbers 2 and 3).
Fourth, the script applies the axis restrictions as mentioned in the *rotation.py* section.
Fifth, the script executes the **rotation3d** function from *rotation.py* and replaces the original coordinate with the new coordinate. The script also prints information about rotation.  

Note: This setup **does** run the risk of rotating the same molecule twice in one iteration, but given that the probability of the same molecule being chosen is overall less than 50% per "movement" step, not to mention that the rotation is completely random either way, we did not find it a significant problem.

#### Movement
When the rot\_not parameter is equal to 0, the script executes movement, with restrictions. Similarly to the rotation, the script defines how many, and which molecule to rotate. It then executes the movement, writing the moved molecules into a new array, mov\_mol\_new.
Next, the movement function checks whether the new molecules are within arbitrarily defined bounds. The bounds were found manually, and rounded accordingly, to ensure that the movement would not push a molecule too far away from the others - this may result in a lower energy, but is not in reality relevant as a result; obviously the easiest way to decrease the potential energy is to lower the field strength, which decreases quadratically with displacement. The restrictions are executed by the following snippet:
~~~
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
    # Reject change:
    write_now = 0
~~~
It is possible to create an opposite transformation in the **else** section of this script, but it was deemed as unneccesary, as the movement in any case is random. Instead, the change is just not accepted, and the script loops back to the next, random step.

#### Writing to files
The script then writes the file using the **writefile** function as described above, given that the parameter writing is given as 1, not 0. If writing is set to 0 (for example in the case of debugging), no writing of the file will occur.

### Optimizing energy
The optimization, quite simply, is done by the *molec.py* functions as desribed above, as long as the parameter optimize is 1. This gives the possibility to decide whether the random states should be optimized or not.
We also put the amount of steps to 100, as this, we found, might increase the chance of localizing a lower energy state.

# Results
When running the simulation on the provided *w6.xyz* files, a series of interesting phenomena occured.
The first thing we noticed was the fact that the simulation very often provided outrageously low energies. [Figure of energies graph] When we inspected this, we found that in some cases, the minimization caused some molecules to fuse together [figure] and others to push the hydrogen atoms closer to the oxygen atoms. Thus, every time we found a low energy point, we inspected the energies manually, and found that in every case where the energy found was lower than the initial energy, there was something wrong in the function.
We then tried to run the simulation again, and in this case, only use the rotation function (defining rot\_not as always 1). Here, we did *not* observe any issue with fused molecules, but we did not observe that the energy was lower either. [graph] This caused us to believe that the *w6.xyz* file was already optimized, and not generated at random, as previously thought. Hence, the problem of movement would impair us greatly, seeing as


# Conclusions

While the study in general proved inconclusive in optimizing the file, we did indeed find a series of interesting issues that would have been interesting to inspect, correct, or test, if given more time to work on the project:
+ Restricting the movement to ensure that no oxygen molecules were less than 3 Å apart. This would, in turn, ensure that there would be no H-H bonds created (as occured), as well as ensure that there would be no oxygen bonded to three hydrogens (which also occured).
+ Restricting rotation such that the entire molecule could rotate collectively around an axis perpendicular to the plane of the atoms.
+ Defining movement using a Lennard-Jones potential, as done by others [cite PELE]
+ Visualizing and exploring the energy landscape of the simulation (for example basinhopping), as done by others  
+ Optimizing the algorithm to run exclusively on openbabel / pybel based code to ensure both compatibility and the script running faster.

In the end, while these trials proved inconclusive, they do open the door for future projects to look into these different areas individually, as the scope of these collectively, far exceeds the time-constraint of a two-week project.   
