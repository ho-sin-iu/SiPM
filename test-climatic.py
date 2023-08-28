#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/test
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 28, 2023
Version     : 1.0
"""

from instrum import *
from analysis import *

#===============================================================
# SETTINGS
TrialName = "0828"	# an identifier of the trial, e.g. date
dispmode = False	# if True, display the real-time status of the system
testmode = False	# if True, skip the configuration of the serial ports

# Serial ports
OS = "Mac" # {"Linux", "Mac", "Windows"}
Tsens_port = {"Mac": "/dev/tty.usbserial-AM01VS5B", "Linux": "/dev/ttyUSB0", "Windows": "COM0"}
Tctrl_port = {"Mac": "/dev/tty.usbserial-14220", "Linux": "/dev/ttyUSB1", "Windows": "COM1"}

# Temperature --- all values are in degree Celsius.
Tmin = 10.0		# minimum temperature
Tmax = 20.0  	# maximun temperature
Tstep = 2.0     # step in temperature
Tend = 30.0     # temperature to restore
Tsetpts = np.arange(Tmin, Tmax, Tstep) # temperature setpoints
print("Info:    Temperature setpoints:", Tsetpts, "°C")

# Sampling periods
tP		= .2	# sampling period in second, time between two samplings
Pwait	= 12	# time (in unit of tP) required to wait for the temperature to change until it stabilizes
Pstab 	= 20	# time (in unit of tP) required to stablize the temperature during the electrical measurements 

# Calculate the ETA
t_esti = tP * (Pwait + Pstab) * len(Tsetpts) + tP * Pstab
t_esti_q, t_esti_r = int(t_esti // 60), t_esti % 60
print(f"Info:    Estimated to take {t_esti_q:d} min {t_esti_r:.2f} sec")

#===============================================================
# MEASUREMENT

# Open the serial ports
Tsens = Sensor(Tsens_port[OS], 115200, 0.1, disp=dispmode, test=testmode, name="Tsens")
Tctrl = Controller(Tctrl_port[OS], 115200, 0.1, disp=dispmode, test=testmode, name="Tctrl")

try:
	# Stepwise heat-up cycle
	for T in Tsetpts:
		Tctrl.set_temp(T)
		Tsens.wait(6, 7, P=Pwait, tsleep=tP)
		Tsens.test_stable(n=Pstab, x=7, Tthr=0.5, tsleep=tP) # Only use the values from sensor T7

	# Restores to room temperature
	Tctrl.set_temp(Tend)
	Tsens.test_stable(n=Pstab, x=7, Tthr=0.5, tsleep=tP) # Only use the values from sensor T7

except KeyboardInterrupt:
	print(f"Error:   KeyboardInterrupt")

# Close the serial ports
Tsens.close() 
Tctrl.close()

np.save(f"./examples/{TrialName}-tT-test", Tsens.tT) # save the data npy file
print("Success: data saved as " + TrialName + "-tT-test.npy")


#===============================================================
# PLOT 

tT = tTcurve()
tT.plot(Tsens.tT)
tT.addLegend()
tT.save(TrialName + "-tT-test.png")
print("Success: Plot saved as " + TrialName + "-tT-test.png")

