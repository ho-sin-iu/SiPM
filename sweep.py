#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/sweep
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 27, 2023
Version     : 2.2
Description :
	- This is the extension of older sweepplot series
    * Performs voltage sweep by ``smu.sweep`` as many times as you wish
    * Analyzes the I-V curve immediately
"""

from instrum import *
from analysis import *
import numpy as np
import pandas as pd

TrialName = "0827"
smu_address = '169.254.0.2' 

#===============================================================
# Measurement

smu = SMU(smu_address, timeout=120000, name='2657A')
smu.write_tsp("TSPbeep")
smu.write_tsp("TSPsweep")
# smu.beep() # let SMU play the music!
# smu.fullname()
iv = smu.sweep(30, 48.0, 54.0, 0.1, 5e-3, trialname="Troom", repeated=1)
smu.close() # remember to turn off the instrument after use!

#===============================================================
# ANALYSIS

# try:
# 	
# except:
# 	print("Error:   no data registered!")

# Plot
IV = IVcurve() 
IV.addScatter(iv[0], iv[1], color="r", label=r"$T_{room}$")

# Obtain breakdown voltage
IV.calcPeak(order=1, pltMode=1, pltItems=["vlines", "text"]) # "vlines", "text"
# vbr1 = IV.calcVbr1(pltItems=["linfit", "polyfit", "span"]) #"linfit", "polyfit", "intersect", "span"
# vbr2 = IV.calcVbr2()
# vbr3 = IV.calcVbr3()
# 
# # Print result
# 
# print(f"Info:    Breakdown voltage V_br\n" + 
#       " " * 9 + f"Method 1 = {vbr1:5.2f} V\n" + 
#       " " * 9 + f"Method 2 = {vbr2:5.2f} V\n" +
#       " " * 9 + f"Method 3 = {vbr3:5.2f} V"
#       )

# Save plot
IV.addLegend(outright=False)
IV.save(name=TrialName + "-IV-test.png")
print("Success: Plot saved as " + TrialName + "-IV-test.png")
