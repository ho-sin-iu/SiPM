# `instrum` branch
- The project `sipmlab` will combine two previous projects, including
  - `temp_test.py` and 
  - `sweepplot_v1.x.x.py` series.
- The programming paradigm will transferred from a procedural one to an object-oriented one.
  - Functions of a device are rewritten as the methods under a class that characterizes that device. 
## 2023-07-18
### `instrum_v1.3.py`
- A module with several classes that describes the instruments.
  - class "SMU"
  - [New] class "SerialPort" and its subclasses Controller and Sensor
- The function "readfile" can import text file that contains instructions for troubleshooting.
  - `helpIP.txt`
  - [New] `helpSerialPermission.txt`
- Output messages are cateforized into `Error`, `Warning`, `Info`, etc. with unified formats.

## 2023-07-25
### `instrum_v1.4.py`
