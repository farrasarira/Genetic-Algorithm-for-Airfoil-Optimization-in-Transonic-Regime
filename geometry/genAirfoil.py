

import numpy as np

def calcMaxThickness(k):
    [x_, y_upper, y_lower] = parsec(k)
    return y_upper.max() - y_lower.min()

def genGeom(k):
    [x_, y_upper, y_lower] = parsec(k)

    df1 = np.column_stack((x_[::-1], y_upper[::-1]))
    df2 = np.column_stack((x_, y_lower))

    with open(".\geometry\current_geom.dat", "w") as f:
        f.write("current_mesh\n")
    with open(".\geometry\current_geom.dat", "ab") as f:
        np.savetxt(f, df1, fmt='  %5f')
    with open(".\geometry\current_geom.dat", "ab") as f:
        np.savetxt(f, df2, fmt='  %5f')

def parsec(k):
    n_pts = 140
    x_pts = (1 - np.cos(np.linspace(0, 1, int(np.ceil(n_pts/2)))*np.pi)) / 2

    # x0 = [0.0147, 0.3015, 0.0599, -0.4360, 0.2996, -0.060, 0.4406, 0, 0.00, 0, 14.67] #NACA0012
    # x0 = [0.083, 0.4312, 0.0629, -0.4273, 0.3441, -0.0588, 0.7018, 0, 0, -6.86, 8.08] #RAE2882
    # x0 = [0.0095, 0.4468, 0.1194, -0.7844, 0.2774, -0.0341, 0.1810, 0, 0, -18.82, 9.74] #RAE2882
    x0 = [0.0090, 0.430, 0.063, -0.4, 0.350, -0.059, 0.73, 0, 0.0025, -12.1, -3.0]

    foil = Airfoil.coordinate_Z(x0)

    zup_pts = np.zeros(len(x_pts))
    zlow_pts = np.zeros(len(x_pts))
    for j in range(len(x_pts)):
        zup_pts[j] = foil.Z_up(x_pts[j])
        zlow_pts[j] = foil.Z_low(x_pts[j])

    return x_pts, zup_pts, zlow_pts
    
    

class Airfoil:
  def __init__(self, parsec_params, Z_up, Z_low):
    self.parsec_params = parsec_params
    self.Z_up = Z_up
    self.Z_low = Z_low
    self.labels = ['r_le', 'X_up', 'Z_up', 'Z_xxup', 'X_low', 'Z_low',
                   'Z_xxlow', 'Z_te', 'dZ_te', 'alpha_te', 'beta_te']

  @staticmethod
  def coefficients(args):
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 = args

    p10, p11 = np.deg2rad([p10, p11])

    coeff = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]

    # Upper surface
    C_up = np.ones((6, 6), dtype=np.float64)
    for i in range(6):
      C_up[1, i] = p2 ** coeff[i]
      C_up[2, i] = coeff[i]
      C_up[3, i] = coeff[i] * p2 ** (coeff[i] - 1)
    C_up[4, 0] = -0.25 * p2 ** (-1.5)
    C_up[4, 1] = 0.75 * p2 ** (-0.5)
    C_up[4, 2] = 3.75 * p2 ** (0.5)
    C_up[4, 3] = 8.75 * p2 ** (1.5)
    C_up[4, 4] = 15.75 * p2 ** (2.5)
    C_up[4, 5] = 24.75 * p2 ** (3.5)
    C_up[5, 1:] = 0

    B_up = np.zeros(6)
    B_up[0] = p8 + p9 / 2
    B_up[1] = p3
    B_up[2] = np.tan(p10 - p11 / 2)
    B_up[4] = p4
    B_up[5] = np.sqrt(2 * p1)

    A_up = np.linalg.solve(C_up, B_up)

    # Lower surface
    C_low = np.ones((6, 6), dtype=np.float64)
    for i in range(6):
      C_low[1, i] = p5 ** coeff[i]
      C_low[2, i] = coeff[i]
      C_low[3, i] = coeff[i] * p5 ** (coeff[i] - 1)
    C_low[4, 0] = -0.25 * p5 ** (-1.5)
    C_low[4, 1] = 0.75 * p5 ** (-0.5)
    C_low[4, 2] = 3.75 * p5 ** (0.5)
    C_low[4, 3] = 8.75 * p5 ** (1.5)
    C_low[4, 4] = 15.75 * p5 ** (2.5)
    C_low[4, 5] = 24.75 * p5 ** (3.5)
    C_low[5, 1:] = 0

    B_low = np.zeros(6)
    B_low[0] = p8 - p9 / 2
    B_low[1] = p6
    B_low[2] = np.tan(p10 + p11 / 2)
    B_low[4] = p7
    B_low[5] = -np.sqrt(2 * p1)

    A_low = np.linalg.solve(C_low, B_low)

    return A_up, A_low

  @classmethod
  def coordinate_Z(cls, *args):
    A_up, A_low = Airfoil.coefficients(*args)

    def Z_up(x):
      _Z_up = 0
      for i in range(6):
        _Z_up += A_up[i] * x ** (i + 0.5)
      return _Z_up

    def Z_low(x):
      _Z_low = 0
      for i in range(6):
        _Z_low += A_low[i] * x ** (i + 0.5)
      return _Z_low

    return Airfoil(*args, Z_up, Z_low)