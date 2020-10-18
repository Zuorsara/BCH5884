#!/usr/bin/env python3
#https://github.com/Zuorsara/BCH5884
#Justin Lee's Project 1

import sys
import math

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

pdbfilename=sys.argv[1]

f=open(pdbfilename)
lines=f.readlines()
f.close()

records=[]
massdict={"H":1.01, "C":12.01,"N":14.01,"O":16.0, "P":30.97, "S":32.07,"MG":24.30}

for line in lines:
    Atom=str(line[0:4])
    Serial=str(line[7:11])
    Name=str(line[13:16]).strip()
    Residue=str(line[17:20])
    Character=str(line[21])
    Sequence=int(line[22:26])
    Xcoord=float(line[30:38])
    Ycoord=float(line[39:47])
    Zcoord=float(line[47:54])
    Occupancy=float(line[55:60])
    Temperature=float(line[61:66])
    Element=line[76:78].strip()
    Mass=massdict[Element]
    records.append([Atom,Serial,Name,Residue,Character,Sequence,Xcoord,Ycoord,Zcoord,Occupancy,Temperature,Element,Mass])

#print (records)



#Placeholders so they are defined for the "record in records" calculation
MassTotal=0
XCM=0
YCM=0
ZCM=0


# ?CM denotes the summation of X, Y, or Z coordinate multipled by it's mass
for record in records:
    #print (record[6],record[7],record[8])
    #print (record[12])

    MassTotal+=record[12]
    XCM+=record[6]*record[12]
    YCM+=record[7]*record[12]
    ZCM+=record[8]*record[12]


#print ("These are my mass totals & position totals", (MassTotal),(XCM),(YCM),(ZCM))



#CMC? denotes the center of mass coordinate for X, Y, or Z
CMCX=XCM/MassTotal
CMCY=YCM/MassTotal
CMCZ=ZCM/MassTotal

#print ("There are my center of mass coordinates", (CMCX),(CMCY),(CMCZ))

#Center of Mass
XNC=0
YNC=0
ZNC=0

f=open("output.pdb",'w')
for serial in records:
    output="{0:5}{1:6}{2:5}{3:5}{4:1}{5:3}{6:12.3f}{7:8.3f}{8:8.3f}{9:8.2f}{10:6.3f}{11:>12} {12:6.2f}\n"
    f.write(output.format(serial[0],serial[1],serial[2],serial[3],serial[4],serial[5],serial[6]-CMCX,serial[7]-CMCY,serial[8]-CMCZ,serial[9],serial[10],serial[11],serial[12]))
f.close()

print("Done")
