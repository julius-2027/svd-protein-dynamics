import numpy as np
from matplotlib import pyplot as plt

U = np.loadtxt('U.txt')
S = np.loadtxt('S.txt')
V_T = np.loadtxt('V_T.txt')

plt.rcParams.update({'font.size': 17})
plt.subplots_adjust(left=0.2,bottom=0.2,top=0.8)
plt.grid(True)

# percent of variance explained by including the first k singular values
plt.title('Percent of variation in protein structure \n explained by rank k approximation')
plt.xlabel('k')
plt.ylabel('100% $\cdot$ $(\sum_{i \leq k} {\sigma_i}^2)$ $\cdot$ $(\sum_{i} {\sigma_i}^2)^{-1}$')
plt.xticks([1,200,400,600,800,S.size])
percent_var = 100*np.cumsum(S**2)/np.sum(S**2)

for k in range(1,10+1):
    print(k, '&', f'{percent_var[k-1]:.2f} \%', '\\')

plt.plot(np.arange(1,S.size+1),percent_var)
plt.savefig('percent_var.pdf')
plt.savefig('percent_var.png')

plt.show()
