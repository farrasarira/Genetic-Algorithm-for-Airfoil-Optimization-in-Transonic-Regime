import argparse
import math
import sys
from pathlib import Path
import subprocess
import os
import time
import fileinput
import shutil

def generateMesh():
    # delete mesh
    # deleteMesh_command = (["rm","current_mesh.su2"])
    # subprocess.run(deleteMesh_command,cwd="./genMesh")

    # generate mesh file (.obj) from geometry data (.dat)
    genMesh_command = ([".\genMesh\mesh.exe", "0.05" , "0.000001" , "1.3" , "2.0" , ".\geometry\current_geom.dat"])
    subprocess.run(genMesh_command)
    
    # convert mesh format from (.geo) to (.su2)
    convMesh_command =([".\genMesh\gmsh.exe",".\current_mesh.geo"])
    p = subprocess.Popen(convMesh_command)
    time.sleep(2.5)
    p.terminate()

    # delete Fluid marker
    for line in fileinput.input(".\current_mesh.su2", inplace=1):
        if "NMARK" in line:
            line = line.replace("NMARK= 3", "NMARK= 2")
        sys.stdout.write(line)    

    #  copy mesh file
    moveMesh_command = (["xcopy",".\current_mesh.su2" , ".\su2","/Y"])
    subprocess.run(moveMesh_command)

    # deleteMesh_command = (["rm","current_mesh.geom"])
    # subprocess.run(deleteMesh_command,cwd="./genMesh")



