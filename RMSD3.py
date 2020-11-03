#!/usr/bin/env python3
#https://github.com/Zuorsara/BCH5884
#Justin Lee's RMSD Homework

import sys
import math
import numpy as np

#Line Formatting
#Colum - Typ - Description
#01-06 - str - ATOM
#07-11 - str - Serial #
#13-16 - str - Name
#17-17 - str - Locator [Didn't Add, Empty Space]
#18-20 - str - Residue Name
#22-22 - str - Character
#23-26 - int - Residue Sequence #
#27-27 - str - iCode? [Didn't Add, Empty Space]
#31-38 - flo - X coordinate
#39-46 - flo - Y coordinate
#47-54 - flo - Z coordinate
#55-60 - flo - Occupancy? [Always seems to be 1.00]
#61-66 - flo - Temperature Factor
#77-78 - flo - Elemental Symbol
#79-80 - ??? - Atomic Charge [Didn't Add, No Charges]]
#Original line #'s were changed to remove the need to strip each one.

#v=sys.argv[1]
#w=sys.argv[2]

def readpdb(pdbfile):
    f=open(pdbfile,'r')
    lines=f.readlines()
    f.close()
    Xcoord=[]
    Ycoord=[]
    Zcoord=[]
    for line in lines:
        words=line.split()
        if words[0]=="ATOM":
            Xcoord.append(float(words[6]))
            Ycoord.append(float(words[7]))
            Zcoord.append(float(words[8]))
    return (Xcoord,Ycoord,Zcoord)

v=readpdb("2FA9noend.pdb")
w=readpdb("2FA9noend2mov.pdb")

def rmsd(v,w):
    Vx=np.array(v[0])
    Vy=np.array(v[1])
    Vz=np.array(v[2])
    Wx=np.array(w[0])
    Wy=np.array(w[1])
    Wz=np.array(w[2])
    Diff1=((Vx-Wx) * (Vx-Wx))
    Diff2=((Vy-Wy) * (Vy-Wy))
    Diff3=((Vz-Wz) * (Vz-Wz))
    n=len(Vx)
    XYZ=[]
    XYZ.append((Diff1)+(Diff2)+(Diff3))
    Answer=math.sqrt(sum(XYZ)/n)

rmsd(v,w)
print(Answer)
