import openbabel as ob
import pybel as pb
import numpy as np
import filecmp
from math import pi, sin, cos

import extract as ex
import rotation as rot
import movement as mov

a, b, c = ex.readfile("w6.xyz") # # reading file

# print "molecule is called", a
# print 'the molecule has ', len(b), 'elements, in the order \n', b
# print 'the coordinates of the corresponding elements are: \n', c
"""
# Test for distances > paste this after rotation to see if the rotation  maintains accurate distance

dlist = []

d = 0
for i in range(3):
    d += c[1][i]-c[0][i]
d = np.sqrt(d)
dlist.append(d)

d = 0
for i in range(3):
    d += c[1][i]-c[2][i]
d = np.sqrt(d)
dlist.append(d)

d = 0
for i in range(3):
    d += c[0][i]-c[2][i]
d = np.sqrt(d)
dlist.append(d)

print "DISTANCE BEFORE: \n", dlist
"""

# # Test to rotate first hydrogen of the molecule around axis of first oxygen (O[1]) and second H
# c[1] = rot.rotation3d(c[0],c[2],c[1], 1*pi)

# # Test to move first OHH set with randommove function
# print b
"""
ohh_1 = c[:3]
# print ohh_1
ohh_1_new = mov.randommove(ohh_1, 5)
# print ohh_1_new

c[:3] = ohh_1_new
"""
c[:3] = mov.randommove(c[:3],5)

# print c

q = ex.divide(b,c)
# print q, '\n', len(q), #'\n', type(q)

c_new = ex.unite(q)
# print c_new

# Test for equality
test = (c==c_new).all()
print test

"""
dlist = []
d = 0
for i in range(3):
    d += c[1][i]-c[0][i]
d = np.sqrt(d)
dlist.append(d)

d = 0
for i in range(3):
    d += c[1][i]-c[2][i]
d = np.sqrt(d)
dlist.append(d)

d = 0
for i in range(3):
    d += c[0][i]-c[2][i]
d = np.sqrt(d)
dlist.append(d)

print "DISTANCE AFTER: \n", dlist
"""

ex.writefile("w6_2.xyz", a, b, c) # # Writing file


"""
# # # Test if identical files
test = filecmp.cmp('w6.xyz','w6_2.xyz')
print test
"""
