from genMesh.genMesh import generateMesh
from geometry.genAirfoil import (genGeom, calcMaxThickness)

import subprocess
import pandas as pd
import numpy as np
import sys
import fileinput

def evaluateAero(indi_denorm,ncon,fun_name):
    nvar = indi_denorm.shape[0] # number of design variables

    if fun_name == "transonic_airfoil":
        # generate airfoil coordinate
        genGeom(indi_denorm[0:-1])

        # generate mesh
        generateMesh()

        # run CFD
        for line in fileinput.input(".\su2\SU2_config.cfg", inplace=1):
            sline=line.strip().split("=")
            if sline[0].startswith("AOA"):
                sline[1]=str(indi_denorm[-1])
            line='='.join(sline)
            print(line)  

        runSU2_command = (["mpiexec","-n","7",".\SU2_CFD.exe",".\SU2_config.cfg"])
        subprocess.run(runSU2_command,cwd=".\su2")

        # get Output
        output = pd.read_csv('.\su2\history.csv')

        cl = output['       "CL"       '].iloc[-1]
        cd = output['       "CD"       '].iloc[-1]
        cmy = output['       "CMz"      '].iloc[-1]
        thickness = calcMaxThickness(indi_denorm[0:-1])

        f = cd      # objective function

        cl_min = 0.824
        cmy_max = 0.093
        thickness_min = 0.12

        g = np.zeros(ncon)
        g[0] = cl_min - cl
        g[1] = abs(cmy) - cmy_max
        g[2] = thickness_min - thickness

    return f, g
    # return f
    

