#!/usr/bin/env python
#https://github.com/Zuorsara/BCH5884.git

import math
import sys

#Module 1: Function that will read in the pdb, parse the information, and separate the information into arrays
	
def readpdb(pdbfilename):
	pdbfile=open(pdbfilename,'r')
	lines=pdbfile.readlines()
	pdbfile.close()
	atom=[]
        number=[]
        segment=[]
        residue=[]
        xcoordinate=[]
        ycoordinate=[]
        zcoordinate=[]
        tempfactor=[]
        element=[]
	for line in lines:
            if lines[0:6]=="ATOM":
		words=line.split()
                atom.append(float(words[0]))
                number.append(float(words[1]))
                segment.append(float(words[2]))
                residue.append(float(words[3]))
                xcoordinate.append(float(words[4]))
                ycoordinate.append(float(words[5]))
                zcoordinate.append(float(words[6]))
                tempfactor.append(float(words[7]))
                element.append(float(words[8]))
            else:
                continue
        return (atom,number,segment,residue,xcoordinate,ycoordinate,zcoordinate,tempfactor,element)

			

def writepdbline(outfile, d):
	"Write a line to outfile in pdb format. d must be a dictionary containing records for an atom"
	outfile.write("%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %1s%-2s\n" % (d['rtype'],d['atomnumber'],d['atomtype'],d['altloc'],d['residue'],d['chain'],d['residuenumber'],d['icode'],d['x'],d['y'],d['z'],d['occupancy'],d['tempfact'],d['element'],d['charge']))

#def findcenter(atomlist):
#	"""Returns the center of a structure. Must specify "geom" for geometric center or "com" for center of mass."""
#	cenx=0
#	ceny=0
#	cenz=0
#	natoms=len(atomlist)
#	for atom in atomlist:
#		cenx+=atom['x']
#		ceny+=atom['y']
#		cenz+=atom['z']
#	cenx=cenx/natoms
#	ceny=ceny/natoms
#	cenz=cenz/natoms
#	return {'x':cenx,'y':ceny,'z':cenz}

def rmsd(atomdict1, atomdict2):
	"""Calculate rmsd between two structures. Structures must have residues in the same order"""
	natoms=len(atomdict1)
	if natoms==len(atomdict2):
		sumsq=0
		for n in range(0,natoms):
			dx=atomdict1[n]['x']-atomdict2[n]['x']
			dy=atomdict1[n]['y']-atomdict2[n]['y']
			dz=atomdict1[n]['z']-atomdict2[n]['z']
			sumsq+=dx*dx+dy*dy+dz*dz
		sumsq=sumsq/natoms
		rmsd=math.sqrt(sumsq)
	else:
		print ("Warning, uneven numbers of atoms. RMSD cannot be calculated")
		rmsd=None
	return rmsd

def getmass(element):
	"""Determine mass for element type"""
	massdict={"H":1.01, "C":12.01,"N":14.01,"O":16.0,"P":30.97,"S":32.07,"MG":24.30}
	mass=massdict.get(element)
	return mass

def findcenter(atomlist, centype="geom"):
	"""Returns the center of a structure. Must specify "geom" for geometric center or "com" for center of mass."""
	cenx=0
	ceny=0
	cenz=0
	natoms=len(atomlist)
	if centype=="geom":
		for atom in atomlist:
			cenx+=atom['x']
			ceny+=atom['y']
			cenz+=atom['z']
		cenx=cenx/natoms
		ceny=ceny/natoms
		cenz=cenz/natoms
	elif centype=="com":
		totmass=0
		for atom in atomlist:
			cenx+=atom['x']*atom['mass']
			ceny+=atom['y']*atom['mass']
			cenz+=atom['z']*atom['mass']
			totmass+=atom['mass']
		cenx=cenx/totmass
		ceny=ceny/totmass
		cenz=cenz/totmass
	else:
		print ("Center type", centype, "not defined")
		sys.exit()
	return {'x':cenx,'y':ceny,'z':cenz}

def centercoords(atomlist,centerdict, check=False):
	"""Centers atoms by subtracting centerdict. Modifies atomlist"""
	for atom in atomlist:
		atom['x']=atom['x']-centerdict['x']
		atom['y']=atom['y']-centerdict['y']
		atom['z']=atom['z']-centerdict['z']
	if check:
		com=findcenter(atomlist,"com")
		geom=findcenter(atomlist,"geom")
		print ("New center of mass is:", com['x'],com['y'],com['z'])
		print ("New geometric center is:", geom['x'],geom['y'],geom['z'])
	return


if __name__=="__main__":
	
	atomlist=readpdb("/home.local/bch5887-01/5tnc_move.pdb")
	atomlist2=readpdb("/home.local/bch5887-01/5tnc_move2.pdb")
	print (atomlist[0])
	r=rmsd(atomlist,atomlist2)
	print ("rmsd =", r)
	geom=findcenter(atomlist, centype="com")
	geom2=findcenter(atomlist2, centype="com")
	print ("geom dict", geom, "geom dict2", geom2)
	centercoords(atomlist,geom,check=True)
	centercoords(atomlist2,geom2,check=True)
	
