from optimization.RGA import (soga)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

npop = 1      # number of population
maxgen = 1     # max generation
pcross = 0.9    # p cross 
pmut = 0.1      # p mutation 

nvar = 12    # number of design variables
ncon = 3    # number of constraints
cht_types = 'G-MCR' # G-MCR , SoF
# eps_const = 1e-4    # equality constraint interval
# bounds for design variables 
#             ['r_le', 'X_up', 'Z_up', 'Z_xxup', 'X_low', 'Z_low','Z_xxlow', 'Z_te', 'dZ_te', 'alpha_te', 'beta_te', 'AOA']
lb = np.array([0.0055, 0.3573, 0.0600, -1.0294, 0.3360, -0.071, -0.0486, -0.02, 0, -20.51,    -3, 2.79])
ub = np.array([0.0215, 0.6043, 0.1194, -0.3900, 0.5376, -0.057,  0.8204,  0.02, 0,     -3, 14.74, 2.79])

# Optimization - Trial
result = soga( "transonic_airfoil", nvar, ncon, lb, ub, npop, maxgen, pcross, pmut, cht_type=cht_types)
best_objval = result[0]
best_individual = result[1]
best_conval = result[2]

print('Best Objective Value: ',best_objval[len(best_objval)-1])
print('Best Individual: ',best_individual[len(best_objval)-1])
print('Best Constraint Value: ',best_conval[len(best_objval)-1])

# ========================================================================
# Plotting the convergence performance
fig = plt.figure()
plt.plot(np.arange(0,maxgen),best_objval)
plt.xlabel('Generation')
plt.ylabel('Obtained Minimum Value')
plt.grid()
plt.show()




