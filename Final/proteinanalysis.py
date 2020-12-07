#!/usr/bin/env python
#https://github.com/Zuorsara/BCH5884.git

#I absolutely bit off more than i could chew, and was unable to get the code to function correctly
#Update 1: pdb files don't load correctly?
#Update 2: the namespace wasn't previously correct

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
		if line[:4]=="ATOM":
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

#Module 3: calculates the total mass of the wildtype and mutant pdb files + other variables for later usage
def findtotalmass():
	sequencenum=[]
	aresidues=[]
	aelements=[]
	totalmass=0
	for record in records:
		sequencenum=np.array(record[1])
		aresidues=np.array(record[3])
		aelements=np.array(record[9])
		totalmass+=record[11]
		
	return sequencenum, aresidues, aelements, totalmass
	

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
		print ("Uneven # of atoms listed, range error)
		rmsd=None
		
	return rmsd

#Module 5: Relative Abundance of the elements
def wildtypeelementalabundance():
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
	plt.savefig("Plot1.png",format="png")
	#plt.show()
	
#Module 7: calculating the mean of the temperature factors
def wildtypetempfactor(wildtype):
	wtempfact=[]
	wtempmean=[]
	for record in records:
		wtempfact=np.array(record[9])
		wtempmean=np.mean(wtempfact)
	return wtempfact, wtempmean

def mutanttempfactor(mutant):
	mtempfact=[]
	mtempmean=[]
	for record in records:
		mtempfact=np.array(record[9])
		mtempmean=np.mean(wtempfact)
	return mtempfact, mtempmean

#Module 8: Graphing the temperature factors
def tfgraph():
	plt.plot(sequencenum, wtempfact, color="blue")
	plt.plot(sequencenum, mtempfact, color="red")
	plt.xlabel('Atom Number')
	plt.ylabel('Temperature Factor')
	plt.savefig("Plot2.png",format="png")
	#plt.show()
	
#Module 9: Relative abundance of residues
def wildtyperesidueabundance(wildtype):
	wnumala=0
	wnumarg=0
	wnumasn=0
	wnumasp=0
	wnumasx=0
	wnumcys=0
	wnumglu=0
	wnumgln=0
	wnumglx=0
	wnumgly=0
	wnumhis=0
	wnumile=0
	wnumleu=0
	wnumlys=0
	wnummet=0
	wnumphe=0
	wnumpro=0
	wnumser=0
	wnumthr=0
	wnumtrp=0
	wnumtyr=0
	wnumval=0
	for residues in aresidues:
		if residue=="ALA":
			wnumala+=1
		elif residue=="ARG":
			wnumarg+=1
		elif residue=="ASN":
			wnumasn+=1
		elif residue=="ASP":
			wnumasp+=1
		elif residue=="ASX":
			wnumasx+=1
		elif residue=="CYS":
			wnumcys+=1
		elif residue=="GLU":
			wnumglu+=1
		elif residue=="GLN":
			wnumgln+=1
		elif residue=="GLX":
			wnumglx+=1
		elif residue=="GLY":
			wnumgly+=1
		elif residue=="HIS":
			wnumhis+=1
		elif residue=="ILE":
			wnumile+=1
		elif residue=="LEU":
			wnumleu+=1
		elif residue=="LYS":
			wnumlys+=1
		elif residue=="MET":
			wnummet+=1
		elif residue=="PHE":
			wnumphe+=1
		elif residue=="PRO":
			wnumpro+=1
		elif residue=="SER":
			wnumser+=1
		elif residue=="THR":
			wnumthr+=1
		elif residue=="TRP":
			wnumtrp+=1
		elif residue=="TYR":
			wnumtyr+=1
		elif residue=="VAL":
			wnumval+=1
		else:
			continue
			
	return wnumala, wnumarg, wnumasn, wnumasp, wnumasx, wnumcys, wnumglu, wnumgln, wnumglx, wnumgly, wnumhis, wnumile, wnumleu, wnumlys, wnummet, wnumphe, wnumpro, wnumser, wnumthr, wnumtrp, wnumtyr, wnumval

def mutantresidueabundance(mutant):
	mnumala=0
	mnumarg=0
	mnumasn=0
	mnumasp=0
	mnumasx=0
	mnumcys=0
	mnumglu=0
	mnumgln=0
	mnumglx=0
	mnumgly=0
	mnumhis=0
	mnumile=0
	mnumleu=0
	mnumlys=0
	mnummet=0
	mnumphe=0
	mnumpro=0
	mnumser=0
	mnumthr=0
	mnumtrp=0
	mnumtyr=0
	mnumval=0
	for residues in aresidues:
		if residue=="ALA":
			mnumala+=1
		elif residue=="ARG":
			mnumarg+=1
		elif residue=="ASN":
			mnumasn+=1
		elif residue=="ASP":
			mnumasp+=1
		elif residue=="ASX":
			mnumasx+=1
		elif residue=="CYS":
			mnumcys+=1
		elif residue=="GLU":
			mnumglu+=1
		elif residue=="GLN":
			mnumgln+=1
		elif residue=="GLX":
			mnumglx+=1
		elif residue=="GLY":
			mnumgly+=1
		elif residue=="HIS":
			mnumhis+=1
		elif residue=="ILE":
			mnumile+=1
		elif residue=="LEU":
			mnumleu+=1
		elif residue=="LYS":
			mnumlys+=1
		elif residue=="MET":
			mnummet+=1
		elif residue=="PHE":
			mnumphe+=1
		elif residue=="PRO":
			mnumpro+=1
		elif residue=="SER":
			mnumser+=1
		elif residue=="THR":
			mnumthr+=1
		elif residue=="TRP":
			mnumtrp+=1
		elif residue=="TYR":
			mnumtyr+=1
		elif residue=="VAL":
			mnumval+=1
		else:
			continue
			
	return mnumala, mnumarg, mnumasn, mnumasp, mnumasx, mnumcys, mnumglu, mnumgln, mnumglx, mnumgly, mnumhis, mnumile, mnumleu, mnumlys, mnummet, mnumphe, mnumpro, mnumser, mnumthr, mnumtrp, mnumtyr, mnumval
		       
#Module 10: Plotting the residue abundance
def residueabundanceplot():
	N = 20
	r1 = (wnumala, mnumala)
	r2 = (wnumarg, mnumarg)
	r3 = (wnumasn, mnumasn)
	r4 = (wnumasp, mnumasp)
	r5 = (wnumcys, mnumcys)
	r6 = (wnumglu, mnumglu)
	r7 = (wnumgln, mnumgln)
	r8 = (wnumgly, mnumgly)
	r9 = (wnumhis, mnumhis)
	r10 = (wnumile, mnumile)
	r11 = (wnumleu, mnumleu)
	r12 = (wnumlys, mnumlys)
	r13 = (wnummet, mnummet)
	r14 = (wnumphe, mnumphe)
	r15 = (wnumpro, mnumpro)
	r16 = (wnumser, mnumser)
	r17 = (wnumthr, mnumthr)
	r18 = (wnumtrp, mnumtrp)
	r19 = (wnumtyr, mnumtyr)
	r20 = (wnumval, mnumval)
	ind = np.arange(N) # the x locations for the groups
	width = 0.075
	fig = plt.figure()
	ax = fig.add_axes([0,0,1,1])
	ax.bar(ind, r1, width, color='r')
	ax.bar(ind, r2, width, color='g')
	ax.bar(ind, r3, width, color='b')
	ax.bar(ind, r4, width, color='r')
	ax.bar(ind, r5, width, color='g')
	ax.bar(ind, r6, width, color='b')
	ax.bar(ind, r7, width, color='r')
	ax.bar(ind, r8, width, color='g')
	ax.bar(ind, r9, width, color='b')
	ax.bar(ind, r10, width, color='r')
	ax.bar(ind, r11, width, color='g')
	ax.bar(ind, r12, width, color='b')
	ax.bar(ind, r13, width, color='r')
	ax.bar(ind, r14, width, color='g')
	ax.bar(ind, r15, width, color='b')
	ax.bar(ind, r16, width, color='r')
	ax.bar(ind, r17, width, color='g')
	ax.bar(ind, r18, width, color='b')
	ax.bar(ind, r19, width, color='r')
	ax.bar(ind, r20, width, color='g')
	ax.set_ylabel('Abundance (in # of elements)')
	ax.set_title('Residue abundance between wildtype and mutant p53 proteins')
	ax.set_xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18', 'G19', 'G20', 'G21', 'G22'))
	ax.set_yticks(np.arange(0, 1500, 100))
	ax.legend(labels=['alanine', 'arginine', 'asparagine', 'aspartic acid', 'cysteine', 'glutamic acid', 'glutamine', 'glycine', 'histidine', 'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine'])
	plt.savefig("Plot3.png",format="png")
	#plt.show()

#Module 11: Opening HTML		       
def openHTML(f,title):
	f.write("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
""")
	f.write("<head>\n")
	f.write("<title>%s</title>\n" % title)
	f.write("</head>\n")
	f.write("<body>\n")

#Module 12: Writing Images		       
def writeHTMLImage(f, title, imgname):
	f.write('<p class="proteininages">%s</p>\n' % title)
	f.write('<img src="%s" />\n' % imgname)
		       
#Module 13: Closing HTML
def closeHTML(f):
	f.write("</body>\n")
	f.write("</html>\n")
	f.close()

#Executing the Code, mutant can be mutants 1, 2, or 3
if __name__=="__main__":
	wildtype=readpdb("'p53 Wildtype.pdb'")
	mutant=readpdb("'p53 Mutant1.pdb'")
	rmsd(wildtype, mutant)
	wildtypeelementalabundance(wildtype)
	mutantelementalabundance(mutant)
	elementalabundanceplot()
	wildtypetempfactor(wildtype)
	mutanttempfactor(mutant)
	tfgraph()
	wildtyperesidueabundance(wildtype)
	mutantresidueabundance(mutant)
	residueabundanceplot()		       
	f=open("output.html",'w')		       
	openHTML(f, "Summary of p53 Protein Wildtype/Mutant Differences")
	f.write("<h1>The average temperature factor of the wildtype p53 is "%s"</h1>" % (wtempmean))
	f.write("<h1>The average temperature factor of the mutant p53 is "%s"</h1>" % (mtempmean))
	f.write("<h1>The RMSD between the wildtype and the mutant is "%s"</h1>" % (rmsd))
	writeHTMLImage(f, "Elemental Abundance", "Plot1.png")
	writeHTMLImage(f, "Temperature Factors", "Plot2.png")
	writeHTMLImage(f, "Residue Abundance", "Plot3.png")	      
	closeHTML(f)
	f.close()
		       
		       
		       
		       
		       
