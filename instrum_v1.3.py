#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/sipmlab
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : July 18, 2023
Version     : v1.3
Description :
	- The project "sipmlab" will combine two previous projects, including
		- temp_test.py and 
		- sweepplot_v1.x.x.py series.
	- The programming paradigm will transferred from a procedural one to an object-oriented one.
		- Functions of a device are rewritten as the methods under a class that characterizes that device. 
	- A module with several classes that describes the instruments.
		- class "SMU"
		* new class "SerialPort" and its subclasses Controller and Sensor
	- The function "readfile" can import text file that contains instructions for troubleshooting.
		- helpIP.txt
		*  helpSerialPermission.txt
	* Output messages are cateforized into "Error", "Warning", "Info", etc. with formats unified.
"""

import numpy as np
#import pandas as pd
import pyvisa
import serial
import time
# import os

#===============================================================

ip_smu = '169.254.0.1' # MacOS
# 192.168.42.234
# 192.168.0.58
# 192.168.32.194
# port_sensor = "/dev/ttyUSB1"
# port_controller = "/dev/ttyUSB0"
port_sensor = "/dev/tty.usbserial-AM01VS5B" # MacOS
port_controller = "/dev/tty.usbserial-14220" # MacOS
disptT = True

#===============================================================

def readfile(fname):
	"""
	Reads a text file (expository texts, data, etc.). 
	Current available files: helpIP.txt, helpSerialPermission.txt

	Parameters
	-------
	fname : string
		file name
	"""
	for line in open(fname):
		print(line, end="")

class SMU():
	def __init__(self, name, ip, timeout): # constructor
		"""

		Activates the SMU via VISA. The correct IP address of the instrument is required.

		Parameters
		--------
		name : string
			SMU's name, e.g. terry. This will be used as the file name when saving the data or plots. 
		ip : string
			IP address of SMU
		timeout : number
			timeout
		"""
		self.name = name
		self.ip = ip
		self.instr = None
		while self.instr == None:
			try:
				rm = pyvisa.ResourceManager()
				rname = 'TCPIP::' + self.ip + '::INSTR' # DON'T USE double quotation mark
				self.instr = rm.open_resource(rname)
				self.instr.timeout = timeout		
				print(f"Success: SMU {self.name} with IP address at {self.ip} installed.")
				break			
			except:
				self.instr = None
				print(f"Error:   The IP address {self.ip} is not found!")
				help = input("Info:    Ask for help? Press 'y' to get possible solutions.")
				if help == "y":
					readfile("helpIP.txt")
				self.ip = input("Info:    Enter another IP address or enter Ctrl+Z to quit.")

	def fullname(self):
		print(self.instr.query("*IDN?"))
	
	def write_tsp(self, fname):
		"""
		Passes TSP scripts to the SMU via Ethernet and let it run the script on its own. 
		Activate the SMU object before applying other methods such as "sweep".

		Parameters
		-------
		fname : string
			file name
		"""
		if self.instr == None: 
			pass
		else:
			for line in open(fname + ".txt"):
				self.instr.write(line)
			self.instr.write(fname + ".save()")
			
	def beep(self):
		if self.instr == None: 
			pass
		else:
			self.instr.write("TSPbeep()")
		
	def sweep(self, T, vstart, vstop, vstep, ilimit, trial, irange="auto"):
		"""vstart, vstop, vstep, ilimit, trial, irange="auto"
		Writes and saves the script that performs a reverse linear voltage sweep into the SMU

		Parameters
		--------
		T : float
			temperature setpoint
		vstart : float
			starting voltage, in volt
		vstop : float
			stopping voltage, in volt
		vstep : float
			step in voltage, in volt
		ilimit : float
			current limit (compliance) of the source, in ampere
		irange : float or string "auto"  
			current measurement range; set either to a string value (in ampere) or to "auto" to enable autorange; defalt: "auto"  
		trial : string
			the number or name of the trial, e.g. "001" or "A1" etc.

		Returns
		-------
		data : string
			data
		"""
		if self.instr == None: 
			pass
		else:
			print("Measuring...")
			self.instr.write("TSPbeep()")
			self.instr.write("TSPsweep()")
			#print(f'sweep({vstart:.2f}, {vstop:.2f}, {vstep:.2f}, {ilimit:.2e}, \"{irange:s}\")')
			self.instr.write(f'sweep({vstart:.2f}, {vstop:.2f}, {vstep:.2f}, {ilimit:.2e}, \"{irange:s}\")')
			#self.instr.write(f'sweep(48.0, 53.0, 0.25, 1e-7, "auto")')		
			data = self.instr.query(f"printData({T:f})")
			print(data)
			# save data?
			return data

	def close(self):
		if self.instr == None: 
			pass
		else:
			self.instr.close()

#===============================================================

class SerialPort():
	def __init__(self, name, port, baudrate, timeout, testmode=False):
		self.name = name
		self.port = port
		self.ser = None
		self.testmode = testmode
		self.starttime = time.time()
		while self.ser == None:
			if self.testmode == True:
				print(f"*** {self.name} in the TEST MODE! ***")
				break
			else:
				try: # configures the port, baud rate, timeout
					self.ser = serial.Serial(self.port)
					self.ser.baudrate = baudrate
					self.ser.timeout = timeout
					print(f"Success: {self.name} with port at {self.port} installed.")			
					break
				except:
					self.ser = None
					print(f"Error:   {self.name} at {self.port} is not found!")
					help = input("Info:    Ask for help? Press 'y' to get possible solutions.")
					if help == "y":
						readfile("helpSerialPermission.txt")
					self.port = input("Info:    Enter another port or enter Ctrl+Z to quit.")

	def get_time(self):
		""" 
		Return the time, which is defined as the time difference between the current time and when the program started.

		Parameters
		--------
		None.

		Returns
		--------
		t : float
			time
		"""
		t = time.time() - self.starttime
		return t

	def close(self):
		if self.testmode == True:
			pass
		else:
			self.ser.close()

class Controller(SerialPort):
    def set_temp(self, T):
        """
        Pass the temperature setpoint to the temperature controller.
		
		Parameters
		--------
		T : float
			temperature in °C

		Returns
		--------
		None.
		"""
        Tbyte = b'set ' + str(T).encode() + b'\r\n'
        #print(Tbyte)
        print("-" * 62 + f"\nSet temperature to T = {T:.2f} °C")
        if self.testmode == True:
            #print("Warning: Temp. controller is not activated!")
            pass
        else:
            self.ser.write(Tbyte) 
            self.ser.read(len(Tbyte))
        return None

class Sensor(SerialPort):
	def get_temp(self, *X):
		"""
		Return the temperature measured by Sensor #6 or #7

		Parameters
		---------
		*X : int or sequence of int
			 sensor number(s)
			 
		Return
		--------
		Tdict : dict
			a dictionary whose keys are the sensor numbers and whose values are the measured temperature values.
		"""
		Tdict = {}
		for x in X: 
			if x not in [6, 7]:
				print("Error:   ValueError---The argument of get_temp must be either 6 or 7.")
				T = np.nan
			else:
				if self.testmode == True: # if the serial port is deactivated 
					#print(f"Warning: Temp. sensor T{x} is not activated!")
					T = np.random.rand() 
				else: 
					self.ser.write(b'tempf? ' + str(x).encode() + b'\r')
					T = float(self.ser.read(7)[0:-2])
			Tdict[x] = T
		#print(Tdict) # uncomment this line to display what is returned
		return Tdict

	def get_tT(self, *X, display_tT=disptT):
		""" 
		Obtain time and temperature. Return the both in a tuple. 
		If the display_tT is toggled on, then it also prints the results in formatted way.

		Parameters
		--------
		 *X : int or sequence of int
			 sensor number(s)

		 display_tT : bool
			whether to print the results
		 
		Return
		-------
		(t, *Tli) : tuple
			a tuple that contains timestamp and temperature values
		"""
		t = self.get_time()
		Tdict = self.get_temp(*X)
		Tli = [Tdict[x] for x in X]
		if display_tT == True:
			text = f"    t = {t:6.2f} sec;      "
			for x in X:
				text += f"T{x} = {Tdict[x]:6.2f} °C;      "
			print(text)
		return (t, *Tli) 

	def wait(self, *X, P, tsleep):
		"""
		Wait for P sample periods while monitoring the temperature.
		
		Parameters
		----------
		*X : int or sequence of int
			sensor number(s)
		P : int
			time (in unit of sample period) required to wait for the temperature to change
		
		Returns
		----------
		None.
		"""
		p = 0 # waiting-cycle counter 
		print(f"Waiting for {P * tsleep:.2f} sec to stablize the temperature...")
		while p < P:
			self.get_tT(*X)
			time.sleep(tsleep) # wait for one sample period
			p += 1

	def test_stable(self, tsleep, n=10, x=7, Tthr=0.2, show_result=False):
		"""
		Executes a loop in which n temperature samples are measured and the stability of temperature is tested. 
		The flag "stable" is initially set as "False". 
		As the loop is executed, if the deviations from the mean of all measured values are less than the set threshold, 
		then the value of "stable" is changed to "True" and the function will jump out of the loop.
		
		Parameters
		--------
		 n : int
		     number of temperature samples, e.g. 10
		 x : int
		     sensor number, either 6 or 7, default = 6
		 Tthr : float
		     temperature fluctuation threshold, default = 0.2 (°C)
		 show_result : bool
		     whether to print the calculation results
		     
		Returns
		-------
		Tmean : float
			average temperature
		"""
		stable = False
		while not stable:
		    # measures n temperature samples 
		    Tarray = np.empty(n)
		    print("Checking stability...") #SCA
		    for i in range(n):
		        Tarray[i] = self.get_tT(x)[1]
		        time.sleep(tsleep) # wait for one sample period
		    
		    # calculates the deviations from the mean.
		    Tmean = Tarray.mean()
		    Tdev = np.abs(Tarray - Tmean) # deviation from mean
		    truthlist = Tdev < Tthr # a list with each entry showing if the deviation from mean is within the threshold 
		    if show_result == True:
		        print(Tarray, Tdev, truthlist)
		    
		    # determines if the temperature is stable
		    # Note: NONE of the values is NOT within the threshold. = ALL of the values ARE within the threshold.
		    [stable, prefix] = [True, ""] if (False not in truthlist) else [False, "UN"]
		    print(f"    T{x} (ave,{n}) = {Tmean:6.2f} °C ............ {prefix}STABLE")
		return Tmean

##===============================================================
# uncomment any of these lines to try the methods of a SMU class.

#smu = SMU('2657A', ip_smu, 10000)
#smu.write_tsp("TSPbeep")
#smu.beep()
#smu.write_tsp("TSPsweep")
#smu.sweep()
#smu.fullname()
#smu.close() # remember to turn off the instrument after use!

##===============================================================
# uncomment the block to test on the whole heat-up cycle

# SERIAL PORT SETTINGS
Tsensor     = Sensor("Temperature sensor", port_sensor, 115200, 0.1, testmode=False)
Tcontroller = Controller("Temperature controller", port_controller, 115200, 0.1, testmode=False)

# TEMPERATURE SETTINGS
# All the temperature values are in °C.
Tmin = 10.0     # minimum temperature
Tmax = 20.0  # maximun temperature
Tstep = 2.5     # step in temperature
Tend = 30     # temperature to restore
Tsetpts = np.arange(Tmin, Tmax, Tstep)

# TIME SETTINGS
tP     = 10     # sample period in second, time between two samples (= measure_interval in Fabian's version)
Pwait  = 6      # time (in unit of sample period) required to wait for the temperature to change until it stabilizes (~ stability_interval in Fabian's version)
Pmeas  = 3      # time (in unit of sample period) required to hold the temperature for the measurement (= niveau in Fabian's version)

#print("\nTemperature setpoints:", Tsetpts, "°C")

# STEPWISE HEAT-UP CYCLE
"""
for i, T in enumerate(Tsetpts):
	Tcontroller.set_temp(T)
	Tsensor.wait(6, 7, P=Pwait, tsleep=tP)
	Tsensor.test_stable(n=10, x=7, Tthr=0.5, tsleep=tP)
"""
# Restores to room temperature
Tcontroller.set_temp(Tend)
Tsensor.test_stable(n=50, x=7, Tthr=0.5, tsleep=tP)

Tsensor.close() 
Tcontroller.close()
print("-" * 62 + "\nEnd")
