# SiPM Characterization

## `instrum.py`
### Updates
#### Older versions
- This module has several classes that controls the instruments.
	- class `SMU`
	- class `SerialPort`
		- subclass `Controller`
		- subclass `Sensor`
- The function `readfile` can import text file that contains instructions for troubleshooting.
	- `helpIPAddress.txt`
	- `helpSerialPermission.txt`
- Output messages are cateforized into **Error**, **Warning**, **Info**, **Success**, etc. with formats unified.

#### 2023-08-28 (v1.6)
* New function `iv_data_parser()` to parse the string into `numpy.ndarray`.
* New local variable `command` in `SMU.sweep` and `SMU.beep`, for debugging.
* Change variable name in `SerialPort.__init__` from `name` to `port`.
* Correct typo and accord the prompt display.
    
## `analysis.py`
### Updates
#### Older versions
- This module analyzes and visualizes the temperature dependence on various electrical properties of SiPM.
- The class `IVcurve` contains several methods that relates to the analysis on the I-V curve:
	- `addPlot`, `addScatter`: overlay a plot on the axes.
	- `Dln` : calculate $\text{d}/\text{d}V[\ln(I)]$ or $\text{d}^2/\text{d}V^2[\ln(I)]$.
	- `Dln` : find the peaks on the curve of $\text{d}/\text{d}V[\ln(I)]$ or $\text{d}^2/\text{d}V^2[\ln(I)]$ vs. $V$.
	- `calcRQ` : calculate the quenching resistance $R_\text{Q}$.
	- `calcVbr1`, `calcVbr2`, `calcVbr3` : calculate the breakdown voltage $V_\text{br}$.
- The class `IVsegment` is specially designed for the method `IVcurve.Vbr1`. It segments a pair of data arrays according to some given range and contains two methods of fit the data:
	- `linfit` : fit the segmented data into the best linear model.
	- `polynfit` : fit the segmented data into the best polynomial model.

#### 2023-08-28 (v1.3)
* The function `labeler` deal with the format to be displayed on a legend.
* Rearrange the plot-related classes.
	* class `Curve` (new) : attributes e.g. `fig`, `ax`; methods `addLegend` and `save`
		* subclass `IVCurve` : new method `labl`, new attribute `name`
		* subclass `tTCurve` (new) : temperature variation with time

## `test-smu.py`
A simple test script to check whether the computer can communicate with the SMU.

### Updates
#### 2023-08-27 (v1.0)
* New subproject created under `SiPM/test`.

#### Examples 
``` 
$ python3 smu-test.py
SMU test starts.
Input the IP address of the SMU:  169.254.0.2
You have successfully connected to [169.254.0.2]
Name of the instrument = [Keithley Instruments Inc., Model 2657A, 1422608, 1.1.6]
Press Y to terminate the connection.  y
Test ends.
```

## `test-sweep.py`
### Updates
#### Older versions
* This is the extension of older `sweepplot` series.
* Performs voltage sweep with SMU by the method `SMU.sweep` 
* Analyzes the $I$-$V$ curve

#### 2023-08-28 (v1.0)


## `test-analysis-IV.py`
The test script branched from original `sweep.py` and older `sweepplot.py`
### Updates
#### 2023-08-28 (v1.0)
* New subproject created under `SiPM/test`.

### Examples
Use the example file `./examples/IV-curve-example.npy` to do the analysis.
For the details of Methods 1, 2 & 3, please refer to the documentation of `analysis.py`.
```
$ python3 test-analysis-IV.py
Info:    Breakdown voltage V_br
         Method 1 = 52.16 V
         Method 2 = 52.33 V
         Method 3 = 51.91 V
Success: Plot saved as 0828-IV-test.png
```
