## rotation.py Version 1.0
## Written by Magnus Berg Sletfjerding (eembees)
###########################################################
"""Rotation Functions, inspired by
http://paulbourke.net/geometry/rotate/
Verified as working 01/07/17 by eembees"""

"""Importing modules"""
import numpy as np
from math import pi, sin, cos, sqrt



def rotation3d(axis1, axis2, point, angle):
    """Rotates a point around an arbitrary axis in a 3D space
    INPUTS:
    *axis1 -- 3d (xyz) array of 1st axis point
    *axis2 -- 3d (xyz) array of 2nd axis point
    *point -- 3d (xyz) array of point to be rotated
    *angle -- angle (radians) of rotation
        Positive angles are counter-clockwise.
    """
    #Translate axis to (virtual) origin, define new point ax
    ax = point - axis1
    axis1_neg = [-x for x in axis1]
    axis2_neg = [-x for x in axis2]

    #Initialize virtual rotation point rot
    rot = [0.0, 0.0, 0.0]

    # Axis direction vector (normalized)
    N = map(sum, zip(axis2, axis1_neg)) # axis vector
    sqsum = (sum(sqrt(N[i]**2) for i in range(3)))
    direction = [N[i]/sqsum for i in range(3)]

    # Simplifying 3d rotation matrix factors - cosine, sine, translation factor
    co = cos(angle)
    tr = (1-cos(angle))
    si = sin(angle)

    x = direction[0]
    y = direction[1]
    z = direction[2]

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


    # # Define rot
    rot[0] = d11 * ax[0] + d12 * ax[1] + d13 * ax[2]
    rot[1] = d21 * ax[0] + d22 * ax[1] + d23 * ax[2]
    rot[2] = d31 * ax[0] + d32 * ax[1] + d33 * ax[2]

    # # Define output point as rotation point transformed back
    newpoint = rot + axis1

    return newpoint

if __name__ == '__main__':
    import extract as ex

    a, b, c = ex.readfile("w6.xyz") # # reading file

    # Test for distances > paste this after rotation to see if the rotation  maintains accurate distance
    '''
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

    dist = np.sqrt(sum([i**2 for i in dlist]))

    print "DISTANCE BEFORE: \n", dist
    '''
    c[1] = rotation3d(c[0], c[2], c[1], 1*pi)

    ex.writefile("w6_rotated.xyz", a, b, c) # # Writing file
