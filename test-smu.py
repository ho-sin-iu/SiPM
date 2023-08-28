#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/test
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 27, 2023
Version     : 1.0
"""
import pyvisa

print('SMU test starts.')
rm = pyvisa.ResourceManager()
success = False

while success != True:
	try:
		addr = input('Input the IP address of the SMU:  ')
		x = rm.open_resource('TCPIP::' + addr + '::INSTR')
		success = True
		print(f'You have successfully connected to [{addr}]')
	except:
		print(f'Sorry! Unable to connect to [{addr}]. Try another IP address!')

whoami = x.query('*IDN?')[:-1]
print(f'Name of the instrument = [{whoami}]')

ask_to_close = ''
while True:
	ask_to_close = input('Press Y to terminate the connection.  ')
	if ask_to_close in ['Y','y']:
		x.close()
		break

print('Test ends.')
