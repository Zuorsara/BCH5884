#!/usr/bin/env python
#https://github.com/Zuorsara/BCH5884.git

import math
import sys
import numpy as np
import matplotlib.pyplot as plt

#Module 1: Function reading in the pdb file, parsing the lines, separating the data 
def readpdb(pdbfilename):
	pdbfile=open(pdbfilename,'r')
	lines=pdbfile.readlines()
	pdbfile.close()
	records=[]
	totalmass=0
	for line in lines:
		if line[:4]=="ATOM"
			data=[]
			data['atom']=line[0:6]
			data['number']=int(line[6:11])
			data['atomtype']=line[12:16]
			data['residue']=line[17:20]
			data['residuenumber']=line[24:28]
			data['x']=float(line[30:38])
			data['y']=float(line[38:46])
			data['z']=float(line[46:54])
			data['occupancy']=float(line[54:60])
			data['tempfact']=float(line[60:66])
			data['element']=line[76:78].strip()
			data['mass']=findmass(data['element'])
			records.append(data)
	
	return records

#Module 2: uses a dictionary of element:mass to find the mass of each element
def findmass(element):
	massdict={"H":1.01, "C":12.01,"N":14.01,"O":16.0,"P":30.97,"S":32.07,"MG":24.30}
	mass=massdict.get(element)
	
	return mass

#Module 3: calculates the total mass of the wildtype and mutant pdb files
def findtotalmass():
	totalmass=0
	for record in records:
		totalmass+=record[11]
		
	return totalmass
	

#Module 4: takes two sets of the previously parsed data for calculating of rmsd, if natoms isn't equal absolute value of the difference is used for the calculate at decreased accuracy
def rmsd(wildtype, mutant):
	nwild=len(wildtype)
	nmut=len(mutant)
	sumsq=0
	if nwild==nmut:
		for n in range(0,nwild):
			dx=wildtype[n]['x']-mutant[n]['x']
			dy=wildtype[n]['y']-mutant[n]['y']
			dz=wildtype[n]['z']-mutant[n]['z']
			sumsq+=dx*dx+dy*dy+dz*dz
		sumsq=sumsq/natoms
		rmsd=math.sqrt(sumsq)
	elif nwild<nmut
		ldiff=nmut-nwild
		lrange=nmut-ldiff
		for n in range(0,lrange):
			dx=wildtype[n]['x']-mutant[n]['x']
			dy=wildtype[n]['y']-mutant[n]['y']
			dz=wildtype[n]['z']-mutant[n]['z']
			sumsq+=dx*dx+dy*dy+dz*dz
		sumsq=sumsq/natoms
		rmsd=math.sqrt(sumsq)
	elif nwild>nmut
		gdiff=nwild-nmut
		grange=nwild-grange
		for n in range(0,grange):
			dx=wildtype[n]['x']-mutant[n]['x']
			dy=wildtype[n]['y']-mutant[n]['y']
			dz=wildtype[n]['z']-mutant[n]['z']
			sumsq+=dx*dx+dy*dy+dz*dz
		sumsq=sumsq/natoms
		rmsd=math.sqrt(sumsq)
	else:
		print ("Warning, uneven numbers of atoms. RMSD cannot be calculated")
		rmsd=None
		
	return rmsd

#Module 5: Relative Abundance of the elements

def wildtypeelementalabundance():
	for record in records:
		sequencenum=np.array(record[1])
		aelements=np.array(record[11])
		
	return sequencenum, aelements

	wnumnit=0
	wnumcar=0
	wnumoxy=0
	for elements in aelements:
		if element=="N":
			wnumnit+=1
		elif element=="C":
			wnumcar+=1
		elif element=="O":
			wnumoxy+=1
		else:
			continue
			
	return wnumnit,wnumcar,wnumoxy

def mutantelementalabundance():
	for record in records:
		aelements=np.array(record[11])
	return aelements
	mnumnit=0
	mnumcar=0
	mnumoxy=0
	for elements in aelements:
		if element=="N":
			mnumnit+=1
		elif element=="C":
			mnumcar+=1
		elif element=="O":
			mnumoxy+=1
		else:
			continue
			
	return mnumnit,mnumcar,mnumoxy

#Module 6: Plotting the elemental abundance (I originally was going to use plotnine, but the tutorial wasn't very clear so I switched back to matplotlib)
def elementalabundanceplot():
	N = 2
	nitrogencomparison = (wnumnit, mnumnit)
	carboncomparison = (wnumcar, mnumcar)
	oxygencomparison = (wnumoxy, mnumoxy)
	ind = np.arange(N) # the x locations for the groups
	width = 0.35
	fig = plt.figure()
	ax = fig.add_axes([0,0,1,1])
	ax.bar(ind, nitrogencomparison, width, color='r')
	ax.bar(ind, carboncomparison, width, color='g')
	ax.bar(ind, oxygencomparison, width, color='b')
	ax.set_ylabel('Abundance (in # of elements)')
	ax.set_title('Elemental abundance between wildtype and mutant p53 proteins')
	ax.set_xticks(ind, ('G1', 'G2',))
	ax.set_yticks(np.arange(0, 1500, 100))
	ax.legend(labels=['Nitrogen', 'Carbon', 'Oxygen'])
	#plt.show()
	
#Module 7: calculating the mean of the temperature factors
def wildtypetempfactor():
	wtempfact=[]
	wtempmean=[]
	for record in records:
		wtempfact=np.array(record[9])
		wtempmean=np.mean(wtempfact)
	return wtempfact, wtempmean

def mutanttempfactor():
	mtempfact=[]
	mtempmean=[]
	for record in records:
		mtempfact=np.array(record[9])
		mtempmean=np.mean(wtempfact)
	return mtempfact, mtempmean

#Module 8: Graphing the temperature factors
def tfgraph():
	plt.plot(sequencenum, absorbance, color="blue")
	plt.plot(sequencenum, 
	plt.xlabel('Time (min)')
	plt.ylabel('Absorbance (mAU)')
	plt.savefig("Plot.png",format="png")
	#plt.show()
