"""This code is inspired by molecule.py by Jan Jensen, as well as code found at
'https://stackoverflow.com/questions/17763655/rotation-of-a-point-in-3d-about-an-arbitrary-axis-using-python'
Importing modules"""
import openbabel as ob
import pybel as pb
import numpy as np
from math import pi, sin, cos

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
        rotated.append(round(sum([r[j][i] * pointToRotate[j] for j in range(3)]),5)) # MAKE SURE TO SPECIFY DECIMAL

    return rotated
