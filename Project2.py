#!/usr/bin/env python3
#https://github.com/Zuorsara/BCH5884.git

import sys
import numpy as np
from scipy.signal import find_peaks
from matplotlib import pyplot as plt



#[1.A] Write a program that will parse the given chromatogram file

data=sys.argv[1]
f=open(data)
lines=f.readlines()
f.close()

time=[]
absorbance=[]
peaks=0
for line in lines[3:]:
    words=line.split()
    try:
        time.append(float(words[0]))
        absorbance.append(float(words[1]))
    except:
        print ("Parsing Complete")
        continue



#[2] Identify the peaks and delineate their boundaries

time=np.array(time)
absorbance=np.array(absorbance)
na=len(absorbance)
peaks, maximum =find_peaks(absorbance, height=100, threshold=None, distance=10)
npeaks=len(peaks)



#[3] For each peak you identify, report the time at which their maximum occured and the maximum absorbance values.

for n in peaks:
    print ("There is a peak at time",time[n],"minutes, at a maximum absorbance value of",absorbance[n],"mAU")



#[1.B] Write a program that will plot the given chromatogram file

plt.plot(time, absorbance)
plt.xlabel('Time (min)')
plt.ylabel('Absorbance (mAU)')
#plt.show()
#plt.savefig("Plot.png",format="png")
