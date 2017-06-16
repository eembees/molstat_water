## plotting.py Version 1.0
## By Magnus Berg Sletfjerding
##############################################

'''
Importing modules
'''

import matplotlib.pyplot as plt
from math import isnan

'''
Importing file data
'''

en_list = []

e = open('energy_list.txt', 'r')

for line in e:
    en_list.append(float(line))

ene_list = [x for x in en_list if isnan(x) == 0]


hist = 1
'''
Plotting
'''
if hist == 1:
    plt.hist(ene_list)
    plt.title("Histogram of MMFF94 Energies")
    filename = 'energy_list_hist.png'
    plt.xlabel("Energy value")
    plt.ylabel("Frequency")
else:
    plt.plot(ene_list)
    plt.title("Graph of MMFF94 Energies")
    filename = 'energy_list_graph.png'
    plt.ylabel("Energy value")
    plt.xlabel("Step")

# plt.show()
plt.savefig(filename)
