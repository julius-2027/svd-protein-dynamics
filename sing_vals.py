import numpy as np
from matplotlib import pyplot as plt

U = np.loadtxt('U.txt')
S = np.loadtxt('S.txt')
V_T = np.loadtxt('V_T.txt')

plt.rcParams.update({'font.size': 20})
plt.subplots_adjust(left=0.2,bottom=0.2)
plt.grid(True)

# singular values
plt.title('Singular values of A')
plt.xlabel('i')
plt.ylabel('$\sigma_i$')
plt.xticks([1,200,400,600,800,S.size])
plt.plot(np.arange(1,S.size+1),S)

plt.savefig('sing_vals.pdf')
plt.savefig('sing_vals.png')

plt.show()
