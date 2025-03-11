import numpy as np
from matplotlib import pyplot as plt
import MDAnalysis as mda
import warnings
warnings.filterwarnings('ignore')

u = mda.Universe('DESRES-Trajectory_sarscov2-11730054-no-water-no-ion-align-CA.psf',
                 'sarscov2-11730054-no-water-no-ion-align-CA-0000.dcd',
                 'sarscov2-11730054-no-water-no-ion-align-CA-0001.dcd',
                 'sarscov2-11730054-no-water-no-ion-align-CA-0002.dcd',
                 'sarscov2-11730054-no-water-no-ion-align-CA-0003.dcd',
                 'sarscov2-11730054-no-water-no-ion-align-CA-0004.dcd')

sel = u.select_atoms('name CA')

print('Building matrix from simulation trajectory:')

X_T = []
for ts in u.trajectory:
    if (ts.frame % 1000 == 0 or ts.frame == len(u.trajectory)-1):
        print('Frame: ', ts.frame, '   ', 'Time: ', u.trajectory.time,'ps')
    X_T.append(sel.positions.flatten())

X_T = np.array(X_T)
X = X_T.T

# subtract off the mean of each row of X
x_mean = np.mean(X,axis=1) # average structure of protein
x_mean_matrix = np.outer(x_mean,np.ones(len(u.trajectory)))
A = X - x_mean_matrix

print('Done building matrix of shape', A.shape,'. Calculating SVD...')

U, S, V_T = np.linalg.svd(A, full_matrices=False, compute_uv=True)

print('Done calculating SVD.')

print('Writing to output files...')

np.savetxt('U.txt',U)
np.savetxt('S.txt',S)
np.savetxt('V_T.txt',V_T)

sel.write('sel.pdb')
u_sel = mda.Universe('sel.pdb')

# exact trajectory of alpha carbons
with mda.Writer('sel.dcd', u_sel.atoms.n_atoms) as w:
    for row in X_T:
        row = np.reshape(row,(-1,3))
        u_sel.atoms.positions = row
        w.write(u_sel.atoms)

# rank 1 approximation
A_1 = S[0] * np.outer(U[:,0],V_T[0])
X_1 = x_mean_matrix + A_1 # add back in the average structure of the protein
X_1_T = X_1.T
with mda.Writer('sel_apprx.dcd', u_sel.atoms.n_atoms) as w:
    for row in X_1_T:
        row = np.reshape(row,(-1,3))
        u_sel.atoms.positions = row
        w.write(u_sel.atoms)

print('Done writing to output files.')
