#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/tempcontrol
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 27, 2023
Version     : 1.0
Description :
	- 
"""

from instrum import *
from analysis import *
import numpy as np
import matplotlib.pyplot as plt

TrialName = "0827"
port_sensor = "/dev/tty.usbserial-AM01VS5B" # MacOS
port_controller = "/dev/tty.usbserial-14220" # MacOS
disp = True

# SERIAL PORT SETTINGS
Tsensor     = Sensor("Temperature sensor", port_sensor, 115200, 0.1, displaymode=disp, testmode=False)
Tcontroller = Controller("Temperature controller", port_controller, 115200, 0.1,  displaymode=disp, testmode=False)

# TEMPERATURE SETTINGS
# All the temperature values are in °C.
Tmin = 10.0		# minimum temperature
Tmax = 20.0  	# maximun temperature
Tstep = 2.0     # step in temperature
Tend = 30.0     # temperature to restore
Tsetpts = np.arange(Tmin, Tmax, Tstep)

# TIME SETTINGS
tP     = .2     # sample period in second, time between two samples
Pwait  = 12      # time (in unit of sample period) required to wait for the temperature to change until it stabilizes
Pmeas  = 0      # time (in unit of sample period) required to hold the temperature for the measurement
Pstab = 20
Pstabend = 40

print("Info:    Temperature setpoints:", Tsetpts, "°C")
t_esti = tP * (Pwait + Pmeas + Pstab) * len(Tsetpts) + tP * Pstabend
t_esti_q, t_esti_r = int(t_esti // 60), t_esti % 60
print(f"Info:    Estimated to take {t_esti_q:d} min {t_esti_r:.2f} sec")

try:
	# STEPWISE HEAT-UP CYCLE  
	for T in Tsetpts:
		Tcontroller.set_temp(T)
		Tsensor.wait(6, 7, P=Pwait, tsleep=tP)
		#measure Pmeas
		Tsensor.test_stable(n=Pstab, x=7, Tthr=0.5, tsleep=tP)

	# Restores to room temperature
	Tcontroller.set_temp(Tend)
	Tsensor.test_stable(n=Pstabend, x=7, Tthr=0.5, tsleep=tP)

except KeyboardInterrupt:
	print(f"Error:   KeyboardInterrupt")

# print(Tsensor.tT)
tT = tTcurve()
tT.plot(Tsensor.tT)

Tsensor.close() 
Tcontroller.close()

tT.addLegend()
tT.save(TrialName + "-tT-test.png")
print("Success: Plot saved as " + TrialName + "-tT-test.png")

