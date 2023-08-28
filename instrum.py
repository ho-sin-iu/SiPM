#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/instrum
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 28, 2023
Version     : 1.6
"""

import numpy as np
import pyvisa 
import serial
import time

#===============================================================

# Common-used IP address 
# ip_address = 169.254.0.2 # uncomment to test locally
# 192.168.42.234
# 192.168.0.58
# 192.168.32.194

# port_sensor = "/dev/ttyUSB1"
# port_controller = "/dev/ttyUSB0"

#===============================================================

def readfile(fpath):
	"""
	Read a text file (expository texts, data, etc.). 

	Parameters
	-------
	fpath : str
		File path.
        
    Returns
    -------
    None.
    
	"""
	for line in open(fpath):
		print(line, end="")

def iv_data_parser(txt):
	iv = np.empty((2,0))
	for i, line in enumerate(txt.splitlines()[1:-1]):
		line_split = line.split('\t') # three strings per line, else: blank line
		if len(line_split) == 3:
			iv_pair = np.array([float(line_split[1]), float(line_split[2])])
			iv = np.column_stack((iv, iv_pair)) # shape = (2,n), col1 = voltage, col2 = current
	return iv

	
#===============================================================

class SMU():
	def __init__(self, address, timeout=10000, name=""):
		"""
		Class constructor. Activate the SMU via VISA.
	
		The attribute ``instr`` is an <fill-in> object.

		Parameters
		----------
		address : str
			IP address of SMU. See ``troubleshooting/helpIPAddress.txtAddress.txt.txt`` for troubleshooting.
		timeout : float, optional
			Timeout. The default is 10000.
		name : str, optional
			SMU's name. The default is ``""``.
		"""
		self.address = address
		self.instr = None
		self.name = name
		while self.instr == None: # Try until success
			try:
				rm = pyvisa.ResourceManager()
				rname = 'TCPIP::'+ self.address +'::INSTR'
				self.instr = rm.open_resource(rname)  # DON'T USE double quotation mark
				self.instr.timeout = timeout		
				print(f"Success: communicating with [SMU {self.name}] at [{self.address}].")
				self.instr.write("smua.nvbuffer1.clear()")
				self.instr.write("smua.nvbuffer2.clear()")
				print(f"Info:    reading buffer cleared.")
				break			
			except:
				self.instr = None
				print(f"Error:   The IP address [{self.address}] is not found!")
				help = input("Info:    Ask for help? Press 'y' to get possible solutions.")
				if help == "y":
					readfile("./troubleshooting/helpIPAddress.txt")
				self.address = input("Info:    Enter another IP address or enter Ctrl+Z to quit.")

	def fullname(self):
		print(self.instr.query("*IDN?"))

	def write_tsp(self, tsp_fname):
		"""
		Pass a TSP script to the SMU via Ethernet, rewrite it and save it in the non-volatile memory of the SMU. 
		Activate the SMU object before applying other methods such as "sweep".

		Parameters
		----------
		fname : string
			file name of the TSP script.
	
		Returns
		-------
		None.
		"""
		if self.instr == None: 
			pass
		else:
			for line in open("./TSP_scripts/" + tsp_fname + ".txt"):
				self.instr.write(line)
			self.instr.write("./TSP_scripts/" + tsp_fname + ".save()")
		print(f"Success: TSP script [{tsp_fname:s}] written.")

	def beep(self, repeated=1, score="7-11"):
		"""
		Run the script called ``TSPbeep`` stored in the the non-volatile memory of the SMU. 
		
		Parameters
		----------
		repeated : int
			number of repeated beeps
		
		Returns
		-------
		None.
		"""
		if self.instr == None: 
			pass
		else:
			self.instr.write("TSPbeep()") # load the script
			command = f'beep({repeated:d}, "{score:s}")'
			# print(command) # debug: do we pass the correct command to the SMU?
			self.instr.write(command)

	def sweep(self, T, vstart, vstop, vstep, ilimit, trialname, irange="auto", repeated=1):
		"""
		Writes and saves the script that performs a reverse linear voltage sweep into the SMU

		Parameters
		----------
		T : float
			Temperature setpoint.
		vstart : float
			Starting voltage, in volt.
		vstop : float
			Stopping voltage, in volt.
		vstep : float
			Step in voltage, in volt.
		ilimit : float
			Current limit (compliance) of the source, in ampere.
		trialname : str
			the number or name of the trial, e.g. "001" or "A1" etc.
		irange : float or string "auto"  
			Current measurement range; set either to a string value (in ampere) or to "auto" to enable autorange; defalt: "auto".  
		repeated : int
			number of repeated sweeps
		
		Returns
		-------
		data : string
			data
			
		Reference
		---------
		.. https://pyvisa.readthedocs.io/en/latest/introduction/event_handling.html#registering-handlers-for-event
		.. https://rfmw.em.keysight.com/wireless/helpfiles/pxivna/Programming/GPIB_Example_Programs/SRQ_Example.htm
		"""
		if self.instr == None: 
			pass
		else:
			try:
				self.beep(1) # beep once at the beginning
				print(f"Info:    voltage sweep from [{vstart:.2f} V] to [{vstop:.2f} V] with step of [{vstep:.2f} V] for [{repeated:d} time(s)]")

				# voltage sweep and request the data
				self.instr.write("TSPsweep()") # load the script
				command = f'sweep({vstart:.2f}, {vstop:.2f}, {vstep:.2f}, {ilimit:.2e}, "{irange:s}", {repeated:d})'
				# print(command) # debug: do we pass the correct command to the SMU?
				self.instr.write(command) # perform voltage sweep
				data = self.instr.query(f"passData({T:f})") # Request the data from reading buffer
				data_np = iv_data_parser(data) # convert the data format from string to numpy.ndarray 
				
				self.beep(2)  # beep twice if completed
				print("Success: measurement completed.")
				
				return data_np
				
			except (KeyboardInterrupt, NotImplementedError):
				print(f"Error:   measurement aborted due to error")
				self.close()
				
	def OPC(self):
		print(self.instr.query("*OPC?"))
		
	def close(self):
		if self.instr == None: 
			pass
		else:
			print(f"Success: [SMU {self.name}] closed.")
			self.instr.close()

#===============================================================

class SerialPort():
	def __init__(self, port, baudrate, timeout, displaymode=False, testmode=False, name=""):
		self.ser = None
		self.displaymode = displaymode
		self.testmode = testmode
		self.name = name
		self.starttime = time.time()
		self.tT = np.empty((3,0))
		while self.ser == None:
			if self.testmode == True:
				print(f"Warning: [{self.name}] in the TEST MODE!")
				break
			else:
				try: # configures the port, baud rate, timeout
					self.ser = serial.Serial(port)
					self.ser.baudrate = baudrate
					self.ser.timeout = timeout
					print(f"Success: [{self.name}] with port at [{port}] installed.")			
					break
				except:
					self.ser = None
					print(f"Error:   [{self.name}] at [{port}] is not found!")
					help = input("Info:    Ask for help? Press 'y' to get possible solutions.")
					if help == "y":
						readfile("troubleshooting/helpSerialPermission.txt")
					port = input("Info:    Enter another port or enter Ctrl+Z to quit.")

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
        if self.displaymode == True:
        	print("-" * 62 + f"\nSet temperature to T = {T:.2f} °C")
        if self.testmode == True:
            print("Warning: the temperature controller is not activated!")
            pass
        else:
            self.ser.write(Tbyte) 
            print(Tbyte)
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
					#print(f"Warning: the temperature sensor T{x} is not activated!")
					T = np.random.rand() 
				else: 
					self.ser.write(b'tempf? ' + str(x).encode() + b'\r')
					T = float(self.ser.read(7)[0:-2])
			Tdict[x] = T
		#print(Tdict) # uncomment this line to display what is returned
		return Tdict

	def get_tT(self, *X):
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
		if self.displaymode == True:
			text = f"    t = {t:6.2f} sec;      "
			for x in X:
				text += f"T{x} = {Tdict[x]:6.2f} °C;      "
			print(text)
		tTli = np.array((t, *Tli))
		#print(tTli)
		self.tT = np.column_stack((self.tT, tTli))
		return tTli

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
		if self.displaymode == True:
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
		idx = {6:1, 7:2}
		while not stable:
			# measures n temperature samples 
			Tarray = np.empty(n)
			if self.displaymode == True:
				print("Checking stability...") #SCA
			for i in range(n): 
				tTli = self.get_tT(6,7)
				Tarray[i] = tTli[idx[x]]
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
			if self.displaymode == True:
				print(f"    T{x} (ave,{n}) = {Tmean:6.2f} °C ............ {prefix}STABLE")
		return Tmean
