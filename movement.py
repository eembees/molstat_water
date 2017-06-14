## movement.py Version 1.0
## Written by Magnus Berg Sletfjerding (eembees)
###########################################################
"""
Movement function:
Moves a set of points along a randomly generated direction vector
INPUTS:
points -- Points as ndarray
dist   -- Distance points move

OUTPUTS:
newpoints -- New Points as ndarray
"""

"""
Importing modules
"""
import numpy as np

"""
Defining move function
"""
def randommove(points, dist):
    # # Standard factors
    dimension = len(points[0])  # # working dimension of coordinate system, extracted from first point
    # # Generate random direction vector
    direction = np.random.random_sample(dimension)
    direction = [(x-0.5)*2 for x in direction]
    # print direction
    sqrtsum_random = np.sqrt(sum([x**2 for x in direction]))
    direction = [direction[i] / sqrtsum_random for i in range(dimension)]

    # # Obtain movement vector
    movement = [direction[i]*dist for i in range(dimension)]

    # # Initiate movement
    newpoints = []              # # List of output points
    # print type(newpoints)

    for i in range(len(points)):
        newpoint = np.add(points[i], movement)
        newpoints.append(newpoint)

    return np.asarray(newpoints)
