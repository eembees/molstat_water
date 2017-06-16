## plotting.py Version 1.0
## By Magnus Berg Sletfjerding
##############################################

'''
Importing modules
'''

import matplotlib.pyplot as plt

'''
Importing file data
'''

ene_list = []

e = open('energy_list.txt', 'r')

for line in e:
    ene_list.append(float(line))

hist = 0
'''
Plotting
'''
if hist == 1:
    plt.hist(ene_list)
    plt.title("Histogram of MMFF94 Energies")
    filename = 'energy_list_hist.png'
else:
    plt.plot(ene_list)
    plt.title("Graph of MMFF94 Energies")
    filename = 'energy_list_graph.png'

plt.xlabel("Energy value")
plt.ylabel("Frequency")
# plt.show()
plt.savefig(filename)
