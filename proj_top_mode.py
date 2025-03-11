import numpy as np
from matplotlib import pyplot as plt

U = np.loadtxt('U.txt')
S = np.loadtxt('S.txt')
V_T = np.loadtxt('V_T.txt')

plt.rcParams.update({'font.size': 20})
plt.subplots_adjust(left=0.2,bottom=0.2,top=0.8)
plt.grid(True)

# movement of protein along the mode of motion given by U[:,0]
plt.title('Projection of protein motion \n onto the top mode')
plt.xlabel('Simulation frame')
plt.plot(V_T[0])

plt.savefig('proj_top_mode.pdf')
plt.savefig('proj_top_mode.png')

plt.show()
