#!/usr/bin/env python3
#https://github.com/Zuorsara/BCH5884
import math

X1 = float(input("Enter X-Coordinate for A: "))
Y1 = float(input("Enter Y-Coordinate for A: "))
X2 = float(input("Enter X-Coordinate for B: "))
Y2 = float(input("Enter Y-Coordinate for B: "))
X3 = float(input("Enter X-Coordinate for c: "))
Y3 = float(input("Enter Y-Coordinate for C: "))

ABx = X1-X2
ABy = Y1-Y2
AB = (ABx*ABx+ABy*ABy)**0.5
BCx = X2-X3
BCy = Y2-Y3
BC = (BCx*BCx+BCy*BCy)**0.5
ACx = X1-X3
ACy = Y1-Y3
AC = (ACx*ACx+ACy*ACy)**0.5

Deg = 57.2958

PreAlpha = ((AB*AB+AC*AC-BC*BC)/(2*AB*AC))
PreBeta = ((BC*BC+AB*AB-AC*AC)/(2*BC*AC))
PreGamma = ((AC*AC+BC*BC-AB*AB)/(2*AC*BC))

Alpha = math.acos(PreAlpha)*Deg
Beta = math.acos(PreBeta)*Deg
Gamma = math.acos(PreGamma)*Deg

print("The Angles of Alpha, Beta, and Gamma are: ",(Alpha),(Beta),(Gamma))
