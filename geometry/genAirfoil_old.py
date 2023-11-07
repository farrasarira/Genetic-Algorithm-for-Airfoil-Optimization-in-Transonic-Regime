"""
Generates PARSEC airfoil.
Get a feel for it: http://www.as.dlr.de/hs/d.PARSEC/Parsec.html
Original paper: http://www.as.dlr.de/hs/h-pdf/H141.pdf
Based on: http://github.com/dqsis/parsec-airfoils
"""

from __future__ import division


from geometry.airfoilgen_baseclass import ParametricAirfoil
import numpy as np


class PARSEC(ParametricAirfoil):

    def __init__(self, k):
        """Takes a dict of coefficients to define PARSEC airfoil.
        Coefficient names: xte, yte, rle, x_suc, y_suc, d2ydx2_suc, th_suc,
                                          x_pre, y_pre, d2ydx2_pre, th_pre"""
        self.k = k        
        try:
            # Parent class contains functions that need to know x-pos of TE
            self.xte = k['xte']
    
            self.coeffs_upper = self._pcoef(k['xte'], k['yte'], k['rle'],
              k['x_suc'], k['y_suc'], k['d2ydx2_suc'], k['th_suc'], 'suction')
            self.coeffs_lower = self._pcoef(k['xte'], k['yte'], k['rle'],
              k['x_pre'], k['y_pre'], k['d2ydx2_pre'], k['th_pre'], 'pressure')
        except TypeError:
            raise Warning("Pass a dict with named coefficients.\n"+
            "Explanation:\n"+self.__init__.__doc__)
        # except KeyError, e:
        #     raise Warning("{:s} was not defined in the dict".format(e))
    
    def __str__(self):
        """Gives some information on airfoil"""
        return ("Airfoil with PARSEC parametrization. Coefficients: {}"
                .format(self.k))

    def _fn_upper_lower(self, xpts):
        return (xpts, self._calc_coords(xpts, self.coeffs_upper),
                xpts, self._calc_coords(xpts, self.coeffs_lower))

    def _calc_coords(self, xpts, coeffs):
        # Powers to raise coefficients to. from __future___ import division!
        pwrs = (1/2, 3/2, 5/2, 7/2, 9/2, 11/2)
        # Make [[1,1,1,1],[2,2,2,2],...] style array
        xptsgrid = np.meshgrid(np.arange(len(pwrs)), xpts)[1]
        # Evaluate points with concise matrix calculations.
        # One x-coordinate is evaluated for every row in xptsgrid
        return np.sum(coeffs*xptsgrid**pwrs, axis=1)

    def _pcoef(self, xte, yte, rle, x_cre, y_cre, d2ydx2_cre, th_cre, surface):
        """Evaluate the PARSEC coefficients.
        From https://github.com/dqsis/parsec-airfoils
        """
        # Initialize coefficients
        coef = np.zeros(6)

        # 1st coefficient depends on surface (pressure or suction)
        if surface == 'pressure':
            coef[0] = -np.sqrt(2*rle)
        elif surface == 'suction':
            coef[0] = np.sqrt(2*rle)

        # Form system of equations
        A = np.array([
                     [xte**1.5, xte**2.5, xte**3.5, xte**4.5, xte**5.5],
                     [x_cre**1.5, x_cre**2.5, x_cre**3.5, x_cre**4.5, 
                      x_cre**5.5],
                     [1.5*np.sqrt(xte), 2.5*xte**1.5, 3.5*xte**2.5, 
                      4.5*xte**3.5, 5.5*xte**4.5],
                     [1.5*np.sqrt(x_cre), 2.5*x_cre**1.5, 3.5*x_cre**2.5, 
                      4.5*x_cre**3.5, 5.5*x_cre**4.5],
                     [0.75*(1/np.sqrt(x_cre)), 3.75*np.sqrt(x_cre), 8.75*x_cre**1.5, 
                      15.75*x_cre**2.5, 24.75*x_cre**3.5]
                     ])
        B = np.array([
                     [yte - coef[0]*np.sqrt(xte)],
                     [y_cre - coef[0]*np.sqrt(x_cre)],
                     [np.tan(th_cre*np.pi/180) - 0.5*coef[0]*(1/np.sqrt(xte))],
                     [-0.5*coef[0]*(1/np.sqrt(x_cre))],
                     [d2ydx2_cre + 0.25*coef[0]*x_cre**(-1.5)]
                     ])
        # Solve system of linear equations
        X = np.linalg.solve(A,B) 
        # Gather all coefficients
        coef[1:6] = X[0:5,0]
        # Return coefficients
        return coef

def genGeom(k):
    # delete airfoil
    # deleteGeom_command = (["rm","current_geom.dat"])
    # subprocess.run(deleteGeom_command,cwd="./geometry")

    
    # Evaluate pressure (lower) surface coefficients
    test_airfoil = PARSEC(k)
    pts = test_airfoil.get_coords()

    df1 = np.column_stack((pts[0][::-1], pts[1][::-1]))
    df2 = np.column_stack((pts[2], pts[3]))

    with open("./geometry/current_geom.dat", "w") as f:
        f.write("current_mesh\n")
    with open("./geometry/current_geom.dat", "ab") as f:
        np.savetxt(f, df1, fmt='  %5f')
    with open("./geometry/current_geom.dat", "ab") as f:
        np.savetxt(f, df2, fmt='  %5f')

    # import matplotlib.pyplot as plt
    # plt.title("PARSEC")
    # plt.plot(pts[0], pts[1], 'o--')
    # plt.plot(pts[2], pts[3], 'o--')
    # plt.gca().axis('equal')
    # plt.show()

def calcMaxThickness(k):
    test_airfoil = PARSEC(k)
    return test_airfoil.max_thickness()


def _example():
    k = {}
    # Sample coefficients
    k['rle'] = 0.01550
    k['x_suc'] = 0.29663
    k['y_suc'] = 0.06002
    k['d2ydx2_suc'] = -0.45150
    k['x_pre'] = 0.29663
    k['y_pre'] = -0.06002
    k['d2ydx2_pre'] = 0.45150
    k['th_suc'] = -8
    k['th_pre'] = 8

    # Trailing edge x and y position
    k['xte'] = 1.0
    k['yte'] = 0.0
    # Evaluate pressure (lower) surface coefficients
    test_airfoil = PARSEC(k)
    pts = test_airfoil.get_coords()

    df1 = np.column_stack((pts[0], pts[1]))
    df2 = np.column_stack((pts[2], pts[3]))

    np.savetxt("current_geom.dat", df1)
    with open("current_geom.dat", "ab") as f:
        np.savetxt(f, df2)
    # f = open("current_geom.dat", "w")
    # f.write("current_geom\n")
    # f.write(pts)

    # import matplotlib.pyplot as plt
    # plt.title("PARSEC")
    # plt.plot(pts[0], pts[1], 'o--')
    # plt.plot(pts[2], pts[3], 'o--')
    # plt.gca().axis('equal')
    # plt.show()


# If this file is run, execute example
if __name__ == "__main__":
    _example()