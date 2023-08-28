#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/test
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 28, 2023
Version     : 1.0
"""

from analysis import *

#===============================================================
# SETTINGS

print_iv = False # if True, print the result before the analysis starts
TrialName = "0828" # an identifier of the trial, e.g. date

#===============================================================
# ANALYSIS

# Import data
ivdata = np.load("./examples/IV-curve-example.npy")
if print_iv == True:
	print(ivdata)

# Plot the I-V curve
IV = IVcurve() 
IV.addScatter(ivdata[0], ivdata[1], color="r")

# Obtain the breakdown voltage
IV.calcPeak(order=1, pltMode=1, pltItems=["vlines", "text"])
vbr1 = IV.calcVbr1(pltItems=["linfit", "polyfit", "intersect", "span"])  
vbr2 = IV.calcVbr2()
vbr3 = IV.calcVbr3()

# Print the result
print(f"Info:    Breakdown voltage V_br")
print(" " * 9 + f"Method 1 = {vbr1:5.2f} V")
print(" " * 9 + f"Method 2 = {vbr2:5.2f} V")
print(" " * 9 + f"Method 3 = {vbr3:5.2f} V")

# Save the plot
IV.addLegend(outright=False)
IV.save(name=TrialName + "-IV-test.png")
print("Success: Plot saved as " + TrialName + "-IV-test.png")
