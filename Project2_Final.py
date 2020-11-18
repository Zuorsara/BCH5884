#!/usr/bin/env python3
#https://github.com/Zuorsara/BCH5884.git

import sys
import numpy as np
from scipy.signal import find_peaks, peak_widths
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
results_half = peak_widths(absorbance, peaks, rel_height = 0.5)
results_adjusted = peak_widths(absorbance, peaks, rel_height = 0.95)
results_full = peak_widths(absorbance, peaks, rel_height = 1.0)
scaling = np.max((time)/na)



#[3] For each peak you identify, report the time at which their maximum occured and the maximum absorbance values.

for i in peaks:
    print ("There is a peak at time",time[i],"minutes, at a maximum absorbance value of",absorbance[i],"mAU")



#[1.B] Write a program that will plot the given chromatogram file

plt.plot(time, absorbance, color="blue")
plt.scatter(results_adjusted[2] * scaling, results_adjusted[1], color="green")
plt.scatter(results_adjusted[3] * scaling, results_adjusted[1], color="green")
plt.xlabel('Time (min)')
plt.ylabel('Absorbance (mAU)')
plt.savefig("Plot.png",format="png")
plt.show()



