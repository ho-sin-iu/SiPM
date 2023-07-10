"""
Author      : Sin-iu Ho (sin-iu.ho@student.uni-tuebingen.de)
Date        : July 10, 2023
Version     : v1.5
Description :
    - remotely control Model 2657A High Power System SourceMeter® Instrument 
	- perform a voltage sweep on a SiPM of Hamamatsu MPPC® S13360 series)
	- measure the voltage and current simultaneously
	- plot the I-V curve
    * calculate the quenching resistance (R_Q) by a linear fit in the forward-biasing region of the I-V curve 
    * evaluate the breakdown voltage (V_BR) using several models
        * Relative derivative method
        * Second derivative method
        * Third derivative method
    * identify the range of Geiger mode

Reference:
    [1] Bonanno et al. (2014). https://doi.org/10.1109/JSEN.2014.2328623 
    [2] Simonetta et al. (2015). https://doi.org/10.1016/j.nima.2015.11.023
    [3] Nagy et al. (2017). https://linkinghub.elsevier.com/retrieve/pii/S0168900217300025
    [4] Hamamatsu MPPC® Technical Note (2022), 
        https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/mppc_kapd9005e.pdf
        page 38: "V_peak is defined as the inflection point of log(I) vs. V curve."
"""

import pyvisa
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from derivative import dxdt
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.optimize import curve_fit

dpi        = 100
measure    = False
version    = "v1.5"

#==============================================================

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

	print("End." + " " * 50)

	keithley.write("smua.source.levelv = 0") 	# set the voltage level back to 0 V
	keithley.write("smua.source.output = 0")	# turn off the source output

	return df

def plotiv(xax, yax, title, ylabel, yscale="linear", plottype="scatter", size=5, color="k", marker="o", lw=2, 
           findpeak=False, linreg=False, legend=False, curvefit=False, ref=""):
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot(111)
    ax.set_xlabel("Voltage $V_{bias}$ [V]")
    ax.set_ylabel(ylabel) # e.g. "Current [A]"
    _title = title + "\n" + ref if ref != "" else title
    ax.set_title(_title)
    plt.yscale(yscale)
    if plottype == "scatter":
        ax.scatter(xax, yax, s=size, c=color, marker=marker, zorder=3, label="data")
    elif plottype == "plot":
        ax.plot(xax, yax, lw=lw, c=color, zorder=3)
    else:
        print("Error: plottype is either \"scatter\" or \"plot\"!")
    
    # for the breakdown voltage measurement V_BR
    #     with Relative Derivative and Second Derivative methods...
    if findpeak == True:
        pk = plotvline(ax, xax, yax, ref)
        if title[0] == "S":
            global Vpeak1, Vpeak2
            Vpeak1 = xax[pk][0]
            #Vpeak2 = xax[pk][1]
   
    # for the breakdown voltage measurement V_BR
    #     with the Third Derivative Method...        
    if curvefit == True:
        xaxx = xax[10:70]
        yaxx = yax[10:70]
        p0 = [3, 0.2, 51, 5] # guess
        global popt, pcov
        popt, pcov = curve_fit(third, xaxx, yaxx, p0, maxfev=5000)
        xrfit = np.linspace(xaxx.min(), xaxx.max(), 100)
        yrfit = third(xrfit, *popt)
        ax.plot(xrfit, yrfit, color="lime", lw=2, zorder=2, label="fitted line")
        pk = plotvline(ax, xrfit, yrfit, ref, h=20)
    
    # for the measurement of quenching resistance R_Q ...
    if linreg == True:
        xs = xax[np.abs(yax) > np.abs(yax.min()) * 0.1]
        ys = yax[np.abs(yax) > np.abs(yax.min()) * 0.1]
        res = linregress(xs, ys)
        xfit = np.linspace(xs.min(), xs.max(), 50)
        yfit = res.intercept + res.slope * xfit
        RQ = 1/res.slope * 1000
        print(f"R_Q    = {RQ:.2f} Ω")
        ax.plot(xfit, yfit, color="m", lw=2, zorder=2, label="fitted line")
        textx = (xs.min() + xs.max())/2 + 0.05
        texty = (ys.min() + ys.max())/2
        text  = f"1/slope = {RQ:.2f} $\Omega$"
        ax.text(textx, texty, text, rotation = 0)
    if legend == True:
        plt.legend()
    
    plt.grid(True, which="both", color="0.9", zorder=1)
    if title[:2] != "$I":
        abbr = title[:2]
    elif "forward" in title:
        abbr = "IVf"
    elif "log" in title:
        abbr = "IVrlog"
    else:
        abbr = "IVr"
    isH = "_HW" if "highwide" in fn_r else ""
    fig.savefig(version + "_" + abbr + isH + ".png")

def plotvline(ax, xax, yax, title, h=2, w=1):
    peaks, _ = find_peaks(yax, height=h, width=w)
    print(f"V_peak = {xax[peaks]} V ... {title:s}")
    for peak in peaks:
        ax.axvline(x = xax[peak], color="k", linestyle="--", zorder=2)
        textx = xax[peak] - 2.5
        texty = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.94
        text  = f"{xax[peak]:.2f} V"
        ax.text(textx, texty, text, rotation = 0)
    return peaks

def third(x, h, sigma, V01, A): # Ref[3]
    coef = A * (2 - h/sigma**2 * (x - V01))
    expon = np.exp(-(x - V01)**2/(2 * sigma**2))
    return coef * expon
    
#==============================================================    

## MAIN
if measure == True:
    ip = '192.168.0.58'
    rm = pyvisa.ResourceManager()
    keithley = rm.open_resource('TCPIP::'+ip+'::INSTR')
    keithley.timeout = 10000 
    
    # forward-biased, find the quenching resistance
    df_f = sweep(0.00, -1.50, -0.025, 1e-1, irange="auto") 
    df_f.to_csv(version + "df_f.csv", index=False)
    
    # reverse-biased, find the breakdown
    df_r = sweep(48.0, 68.0, 0.05, 5e-3, irange="auto")
    df_f.to_csv(version + "df_r.csv", index=False)
    
    # stop sweeping
    keithley.close()
    df_all = pd.concat([df_f, df_r])
    #print(df_all)
else:
    # import existed data
    fn_f = "SiPM_IV_curve_v1.3_forward-biased.csv" #<==============
    df_f = pd.read_csv(fn_f, sep="\t")
    
    fn_r = "SiPM_IV_curve_v1.2_reverse-biased-wide.csv" #<==============
    df_r = pd.read_csv(fn_r, sep="\t")

#==============================================================    
## 1. forward-biased, find the quenching resistance
xf = df_f["Voltage"].to_numpy()
yf = df_f["Current"].to_numpy()

plotiv(xf, yf*1000, "$I$-$V$ curve (forward biasing)", "Current [mA]", size=10, linreg=True, legend=True)
#==============================================================    
## 2. reverse-biased, find the breakdown voltage
xr = df_r["Voltage"].to_numpy()
yr = df_r["Current"].to_numpy()
Dlny  = dxdt(np.log(yr), xr, kind="finite_difference", k=1)
D2lny = dxdt(Dlny, xr, kind="finite_difference", k=1)
D3lny = dxdt(D2lny, xr, kind="finite_difference", k=1)

ys    = [yr, Dlny, D2lny, D3lny]
plottypes = ["scatter", "plot", "plot", "scatter"]
titles  = ["$I$-$V$ curve (reverse biasing, log)", 
           "Relative derivative method",
           "Second derivative method",
           "Third derivative method"]
refs = ["", 
        "(Bonanno et al., 2014)", 
        "(Simonetta et al., 2015)", 
        "(Nagy et al., 2017)" ]
ylabels = ["Current $I$ [A]", 
           "$D\,\ln{(I)}$",
           "$D^2\ln{(I)}$",
           "$D^3\ln{(I)}$"]
yscales = ["log", "linear", "linear", "linear"] 
colors = ["k", "r", "b", "g"]

plotiv(xr, yr*1000, "$I$-$V$ curve (reverse biasing)", "Current [mA]", legend=True)
for i, yri in enumerate(ys):
    find = True if i in [1,2] else False
    legend = True if i in [0,3] else False
    curvefit = True if i in [3] else False
    plotiv(xr, yri, titles[i], ylabels[i], yscales[i], plottype=plottypes[i], 
           color=colors[i], findpeak=find, legend=legend, curvefit=curvefit, size=10,
           ref=refs[i])
print(f"V_BR   = {Vpeak1 - 0.18:.2f} V ... (Hamamatsu MPPC® Technical Note, 2022)") #Ref[4]




