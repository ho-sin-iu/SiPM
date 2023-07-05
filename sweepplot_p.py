"""
Author      : Sin-iu Ho (sin-iu.ho@student.uni-tuebingen.de)
Date        : July 4, 2023
Version     : v1.0
Description :
	* remotely control Model 2657A High Power System SourceMeter® Instrument 
	* perform a voltage sweep on, e.g. a diode or even a Silicon Photomultiplier (SiPM).
	* measure the voltage and current simultaneously and plot the results.
	* suitable for characterization of I-V curve.
"""

import pyvisa
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ip = '192.168.0.58'
rm = pyvisa.ResourceManager()
keithley = rm.open_resource('TCPIP::'+ip+'::INSTR')
keithley.timeout = 10000
#keithley.query('*IDN?')

	
def sweep(vstart, vstop, vstep, ilimit, irange="auto"):
	"""
	DESCRIPTION
		This function performs a reverse linear voltage sweep 

	PARAMETERS
		vstart	: starting voltage, in volts
		vstop	: stopping voltage, in volts
		vstep	: step voltage, in volts
		ilimit	: current limit of the voltage source, in ampere
		irange	: current measurement range; set either to a string value or to "auto" to enable autorange; defalt: "auto"  
	"""

	# reset and initialize
	keithley.write("reset()")
	keithley.write("errorqueue.clear()")
	keithley.write("status.reset()")

	keithley.write("smua.source.func = 1") 			# select voltage as the source function
	keithley.write("display.smua.measure.func = 0") # display the current measurement
	keithley.write("smua.measure.autozero = 1")		# perform autozero once, then disable it
	
	# configure the reading buffer
	keithley.write("smua.nvbuffer1.clear()")
	keithley.write("smua.nvbuffer1.appendmode = 1")
	keithley.write("smua.nvbuffer1.collecttimestamps = 0")
	keithley.write("smua.nvbuffer1.collectsourcevalues = 0")
	keithley.write("smua.nvbuffer1.fillmode = smua.FILL_ONCE")
	keithley.write("smua.nvbuffer2.clear()")
	keithley.write("smua.nvbuffer2.appendmode = 1")
	keithley.write("smua.nvbuffer2.collecttimestamps = 0")
	keithley.write("smua.nvbuffer2.collectsourcevalues = 0")
	keithley.write("smua.nvbuffer2.fillmode = smua.FILL_ONCE")

	keithley.write("smua.measure.autorangev = 1")	# set the voltage measurement range

	# set the current measurement range
	if irange == "auto":
		keithley.write("smua.measure.autorangei = 1")
	else:
		keithley.write("smua.measure.autorangei = 0")
		keithley.write("smua.measure.rangei = " + str(irange))
		
	keithley.write("smua.source.limiti = " + str(ilimit)) # set the current limit (compliance)
	keithley.write("smua.measure.nplc = 1.0") 	# set the integration aperture for measurements to 1 power line cycle (PLC)
	keithley.write("smua.measure.delay = 0.2") 	# set the measurement delay to 0.2 sec

	varray = np.arange(vstart, vstop, vstep)	# arithmetric sequence of voltages to ouput
	df = pd.DataFrame({"Voltage": [], "Current": []}) # dataframe to record data

	print(f"Sweeping from {vstart} V to {vstop} V in steps of {vstep} V...")
	keithley.write("smua.source.output = 1") 	# turn on the source output

	# for-loop that iterates the output values
	for i, vout in enumerate(varray):
		print(f"source voltage = {vout:.3f} V", end="\r")
		keithley.write("smua.source.levelv = " + str(vout)) # update the source level
		keithley.write("smua.measure.iv(smua.nvbuffer1, smua.nvbuffer2)") # measure voltage and current

		# retrieve an print the readings
		ireading = keithley.query("print(smua.nvbuffer1.readings["+ str(i+1) + "])") 
		vreading = keithley.query("print(smua.nvbuffer2.readings["+ str(i+1) + "])")
		ival = float(ireading)
		vval = float(vreading)
		df.loc[len(df)] = [vval, ival]

	#print("end")
	keithley.write("smua.source.levelv = 0") 	# set the voltage level back to 0 V
	keithley.write("smua.source.output = 0")	# turn off the source output

	return df

df1 = sweep(48.00, 50.50, 0.25, 1e-5, irange=10e-9) 
df2 = sweep(50.50, 53.00, 0.10, 1e-5, irange=100e-9)
keithley.close()

print("End." + " " * 50)
df = pd.concat([df1, df2])
print(df)

vax = df["Voltage"].to_numpy()
iax = df["Current"].to_numpy()

fig = plt.figure(dpi=150)
ax = fig.add_subplot(111)
ax.set_xlabel("Voltage [V]")
ax.set_ylabel("Current [A]")
ax.scatter(vax, iax, s=3, c="r", marker="o")
fig.savefig("SiPM_IV_curve.png")



