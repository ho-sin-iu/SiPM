#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/test
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 28, 2023
Version     : 1.0
"""

from instrum import *
import numpy as np

#===============================================================
# OUTPUT SETTINGS
TrialName = "0828"
print_iv = False # if True, print the result when the measurements end.

#===============================================================
# IP ADDRESS SETTINGS
smu_address = '169.254.0.2' 

# An example of manual configuration:
#   IP address of PC    = 169.254.90.209
#   IP address of SMU   = 169.254.0.2
#   Gateway address     = 169.254.255.255
#   Subnet mask         = 255.255.0.0

#===============================================================
# MEASUREMENT

smu = SMU(smu_address, timeout=120000)
smu.write_tsp("TSPbeep")
smu.write_tsp("TSPsweep")
ivdata = smu.sweep(30, 48.0, 64.0, 0.3, 5e-3, trialname="Troom", repeated=1)
smu.close() # ALWAYS turn off the instrument after use!
np.save(TrialName + "-IV-test", ivdata)
