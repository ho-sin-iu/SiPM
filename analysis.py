#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project     : SiPM/analysis
Author      : Sin-iu Ho <sin-iu.ho(at)student.uni-tuebingen.de>
Date        : August 27, 2023
Version     : 1.3
Description :
	- This module analyzes and visualizes the temperature dependence on various electrical properties of SiPM.
    - The class "IVcurve" contains several methods that relates to the analysis on the I-V curve:
        - ``addPlot``, ``addScatter``, ``addLegend`` : overlay a plot on the axes.
        - ``Dln`` : calculate d/dV[ln(I)] or d^2/dV^2[ln(I)].
        - ``Dln`` : find the peaks on the curve of d/dV[ln(I)] or d^2/dV^2[ln(I)] vs. V.
        - ``calcRQ`` : calculate the quenching resistance R_Q.
        - ``calcVbr1``, ``calcVbr2``, ``calcVbr3`` : calculate the breakdown voltage V_br.
    - The class "IVsegment" is specially designed for the method ``IVcurve.Vbr1``. It segments a pair of data arrays according to some given range and contains two methods of fit the data:
        - ``linfit`` : fit the segmented data into the best linear model.
        - ``polynfit`` : fit the segmented data into the best polynomial model.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from derivative import dxdt
from scipy.signal import find_peaks
from scipy.stats import linregress
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

#==============================================================================
class Curve():
    def __init__(self, xlabel, ylabel, title, yscale, grid=True):
        """
        Class constructor. Define the attributes "fig" and "ax" while configuring the basic elements of the plot, such as labels for axes and titles.

        Parameters
        ----------
        xlabel : string, optional
            The label for the x-axis.
        ylabel : string, optional
            The label for the y-axis.
        title : string, optional
            The title for the axes. The default is None.
        yscale : string, optional
            The scale of the y-axis. The default is "linear".
        grid : bool, optional  
            Whether to plot the grid lines.
        
        Returns
        -------
        None.

        """
        self.fig, self.ax = plt.subplots(figsize=(7,4))
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel) 
        self.ax.set_title(title)
        self.ax.set_yscale(yscale)
        plt.grid(grid, which="both", color="0.9", zorder=1)

    def addLegend(self, outright=True):
        """
        Place a legend that is outside or inside the plot

        Parameters
        ----------
        outright : bool, optional
            Whether to place the legend outside. The default is True.

        Returns
        -------
        None.

        """
        if outright == True:
            box = self.ax.get_position()
            self.ax.set_position([box.x0, box.y0, box.width*0.990, box.height])
            plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        else:
            plt.legend(loc="best")
            
    def save(self, name, dpi=200):
    	self.fig.savefig(name, dpi=200)

class IVcurve(Curve):
    def __init__(self, xlabel="Voltage [V]", ylabel="Current [A]", title="$I$-$V$ curve", yscale="log"):
    	Curve.__init__(self, xlabel, ylabel, title, yscale)
	
    def setColors(self, fnli):
        """
        Use the list of file numbers (fnli) to define a list of colors (cli) that creates a color gradient from black (0, 0, 0) to red (1, 0, 0). 
        
        Smaller index in fnli has darker color.

        Parameters
        ----------
        fnli : list of strings
            List of file numbers.

        Returns
        -------
        cli: list
            List of colors in RGB (red, green, blue) tuples of float values.
        """
        cli = [((i+1)/len(fnli), 0.0, 0.0) for i in range(len(fnli))]        
        self.colorset = cli
        return cli
    
    def addPlot(self, vax, iax, color=None, label=""):
        """
        Plot the I-V curve without markers. 
        
        Also pass the values of parameters to define the attributes "vax", "iax", "label", and "color".
        
        Fixed format: 
            * linestyle = "-" (solid line)
            * linewidth = 2 
            * zorder = 3

        Parameters
        ----------
        vax : array-like, shape (n, )
            Voltage data.
        iax : array-like, shape (n, )
            Current data.
        color : array-like or list of colors or color, optional
            The color of marker. The default is None
        label : string, optional
            The label to display in the legend. The default is "".
        
        Returns
        -------
        None.
        """
        self.ax.plot(vax, iax, c=color, lw=2, zorder=3, label=label)
        self.vax = vax
        self.iax = iax
        self.label = label
        self.color = color
        
    def addScatter(self, vax, iax, color=None, label=""):
        """
        Plot the I-V curve as a scatter plot. 
        
        Also pass the values of parameters to define the attributes "vax", "iax", "label", and "color".
        
        Fixed format: 
            * marker = "." (point marker)
            * markersize = 10
            * zorder = 3
            
        Parameters
        ----------
        vax : array-like, shape (n, )
            Voltage data.
        iax : array-like, shape (n, )
            Current data.
        c : array-like or list of colors or color, optional
            The marker colors. Default is None
        label : string, optional
            The label to display in the legend. Default is an empty string.
        
        Returns
        -------
        None.
        """
        self.ax.scatter(vax, iax, color=color, marker=".", s=10, zorder=3, label=label)
        self.vax = vax
        self.iax = iax
        self.label = label
        self.color = color
    
    def Dln(self, order=1, plot=False):
        """
        Evaluate the first or the second derivative of natural logarithm of current with respect to voltage, i.e. d/dV[ln(I)] or d^2/dV^2[ln(I)].
        
        The differentiation is carried out by the functional interface ``derivative.dxdt`` with the method ``FiniteDifference``.

        Parameters
        ----------
        order : int, optional
            Order of the derivative. The default is 1.
        plot : bool, optional
            Whether to plot the derivative vs. V. The default is False.

        Returns
        -------
        D : array-like, shape (n, )
            derivative, either d/dV[ln(I)] or d^2/dV^2[ln(I)]
            
        Raises
        ------
        TypeError 
            If ``order`` is not an integer.
        ValueError
            If ``order`` is neither 1 nor 2.
        
        References
        ----------
        .. [derivative.dxdt] https://derivative.readthedocs.io/en/latest/api.html
        .. [FiniteDifference] https://derivative.readthedocs.io/en/latest/modules.html#finitedifference)

        """
        Dlny = dxdt(np.log(self.iax), self.vax, kind="finite_difference", k=1)
        cdict = {1: "b", 2: "g"}
        sdict = {1: " 1st deri.", 2: " 2st deri."}
        if not isinstance(order, int):
            raise TypeError("order must be an integer.")
        if order not in [1,2]:
            raise ValueError("order can only be 1 (first derivative of ln(I)) or 2 (second derivative of ln(I)).")
        if order==1:
            D = Dlny
        if order==2:
            D = dxdt(Dlny, self.vax, kind="finite_difference", k=1)
        if plot == True:
            self.ax.get_yaxis().set_visible(False)
            self.ax.plot(self.vax, D, c=cdict[order], ls=":", lw=1, zorder=4, label=self.label + sdict[order])
        return D    
    
    def calcPeak(self, order=1, pltMode=0, pltItems=[]):
        """
        Call the method ``IVcurve.Dln`` to evaluate d/dV[ln(I)] or d^2/dV^2[ln(I)] and find the peaks on that derivative by ``scipy.signal.find_peaks``.
        
        The conditions for the peak identification (*EMPIRICAL*) include ``height = 2, width = 1``.

        Parameters
        ----------
        order : int, optional
            Order of the derivative, which is passed to the method "Dln". The default is 1. Acceptable values are
            
                * 1 : find peaks in the 1-d array of d/dV[ln(I)],
                * 2 : find peaks in the 1-d array of d^2/dV^2[ln(I)]
        pltMode : int, optional
            Plotting mode. The default is 0. Acceptable values are
            
                * 0 : act on nothing
                * 1 : act on the first peak
                * 2 : act on all peaks
        pltItems : {"vlines", "text"}, optional
            Items to be displayed on the plot. The default is ``[]``. 
            Following options are supported,
            
                * "vlines" : draw blue dashed vertical lines to indicate peaks
                * "text" : add a text label next to the the vertical line

        Returns
        -------
        first_peak : float or numpy.nan
            The voltage of the leftmost peak (smallest voltage), which is empirically the position of the breakdown voltage and can 
            serve as the breakdown voltage V_br or the initial guess in Method 1 for finding V_br. (see `IVcurve.calcVbr1`)
            
            If no peaks are found, then it returns numpy.nan (Not a Number).
        
        Raises
        ------
        TypeError 
            If ``pltMode`` is not an integer.
        ValueError
            If ``pltMode`` is neither 0 nor 1 nor 2.
            
        References
        ----------
        .. [scipy.signal.find_peaks] https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        
        """
        D = self.Dln(order, plot=False)
        peaks, _ = find_peaks(D, height=2, width=1) # height=2, width=1 --> EMPIRICAL value

        if not isinstance(pltMode, int):
            raise TypeError("pltMode must be an integer.")
        if pltMode not in [0,1,2]:
            raise ValueError("pltMode can only be 0 (act on nothing), 1 (act on the first peak) or 2 (act on all peaks).")
        else:
            modedict = {0: 0, 1: 1, 2: len(peaks)}
            
        for i in range(modedict[pltMode]):
            if "vlines" in pltItems: # draw blue dashed vertical lines to indicate peaks
                self.ax.axvline(x = self.vax[peaks[i]], lw=1, color="b", linestyle="--", zorder=2) 
            if "text" in pltItems: # add a text label next to the the vertical line
                textx = self.vax[peaks[i]] + 0.25
                ybottom, ytop = [self.ax.get_ylim()[j] for j in range(2)]
                texty = ytop * 0.01
                text  = f"{self.vax[peaks[i]]:.2f} V"
                self.ax.text(textx, texty, text, rotation = 90)
        try:
            first_peak = self.vax[peaks][0]
        except:
            first_peak = np.nan
        return first_peak

    def calcRQ(self, N, pltItems=[]):
        """
        Perform the linear regression (on the forward-biasing region) and calculates the quenching resistance R_Q of SiPM.
        
        Parameters
        ----------
        N : int 
            Number of microcells of SiPM.
        pltItems : {"fit", "text"}, optional
            Items to be displayed on the plot. The default is ``[]``. 
            Following options are supported,
            
                * "fit" : plot the linear fit
                * "text" : add a text label next to the the vertical line
  
        Returns
        -------
        RQ : string
            Quenching resistance (in ohm) of the entire SiPM, which is equivalently N quenching resistors in parallel.
        RQerr : string
            Uncertainty of the quenching resistance in ohm.
            

        References
        ----------
        .. [scipy.stats.linregress]  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
        .. [2] https://www.statology.org/standard-error-of-regression-slope/

        """
        vax, iax = self.vax, self.iax
        vslice = vax[np.abs(iax) > np.abs(iax.min()) * 0.1] # slice off the nonlinear part of the I-V curve with current values 0.1 times the maximum --> EMPIRICAL value
        islice = iax[np.abs(iax) > np.abs(iax.min()) * 0.1] # slice off the nonlinear part of the I-V curve with current values 0.1 times the maximum --> EMPIRICAL value
        slope, yintercept, r_value, p_value, std_err = linregress(vslice, islice) #
        vfit = np.linspace(vslice.min(), vslice.max(), len(vslice))
        ifit = yintercept + slope * vfit
        #print(slope, intercept, r_value, p_value, std_err)
        RQ = N * 1/slope #
        RQerr = N * np.sqrt(1/(len(vslice)-2) * ((islice-ifit)**2).sum() / ((vslice - vslice.mean())**2).sum())
        
        if "fit" in pltItems:
            self.ax.plot(vfit, ifit, color="c", lw=3, zorder=2, label=self.label + " (linear fit)")
        if "text" in pltItems:
            textx = (vslice.min() + vslice.max())/2 + 0.05 # right to the middle point of the line. 0.05--> EMPIRICAL value
            texty = (islice.min() + islice.max())/2 - 0.01 # lower to the middle point of the line. 0.01--> EMPIRICAL value
            text  = f"1/slope = {1/slope:.2f} $\Omega$\n$R_Q = N$/slope = {RQ/1e3:.2f} k$\Omega$\n$r$-value = {r_value:.4f}\nstandard error = {std_err:.4f}"
            self.ax.text(textx, texty, text, rotation = 0)
        
        return RQ, RQerr
        
    def calcVbr1(self, color=None, pltItems=[]):
        """
        Calculate the breakdown voltage by finding the intercept of the degree-n fitted polynomial and the "baseline" on semi-log scale. (Intercept method)

        Parameters
        ----------
        color : array-like or list of colors or color, optional
            The color of the plot. The default is None, which will take the same color as the I-V plot.
        pltItems : {"span", "linfit", "polyfit", "intersect"}, optional
            Items to be displayed on the plot. The default is ``[]``. 
            Following options are supported,
            
                * "span" : draw the two vertical spans across the x-axis, representing the baseline (left, cyan) and the working range (right, yellow), respectively.
                * "linfit" : plot the linear fit (without label).
                * "polyfit" : plot the best polynomial fit.
                * "intersect" : plot the intersection of the linear fit and the polynomial fit.

        Returns
        -------
        bestvbr : float
            Breakdown voltage of the best fit.
            
        References
        ----------
        * (Parabolic fitting) \Z. Li *et al.* (2012). *NIMPR, Section A* **695** pp.222-225. https://doi.org/10.1016/j.nima.2011.12.037
        * (Parabolic fitting) N. Dinu *et al.* (2017). *NIMPR, Section A* **845** pp.64-68. https://doi.org/10.1016/j.nima.2016.05.110
        * (Summary of methods) \F. Nagy *et al.* (2017). *NIMPR, Section A* **849** pp.55-59. https://linkinghub.elsevier.com/retrieve/pii/S0168900217300025
        """
        if color == None:
            color = self.color # Take the same color as the I-V curve.
        vax, iax = self.vax, self.iax
        vpeak = self.calcVbr3()  # the initial guess of breakdown voltage given by IVcurve.calcVbr3
        degrange = np.arange(2,7,1) # Try from Degree 2 (quadratic) to Degree 6 polynomial
        
        # Use vpeak to divide the sequence into two segments: 
        BL = IVsegment(vax, iax, vax.min(), vpeak-0.1)  # IVsegment left to vpeak: BL (baseline) ; -0.2 (volts) ---> EMPIRICAL
        WR = IVsegment(vax, iax, vpeak, vax.max()-3)    # IVsegment right to vpeak: WP (working range) ; -3 (volts) ---> EMPIRICAL
        vfit = BL.xfiner() # volage array with the same endpoints as vax and a finer spacing (num = 1000)
        _iVpeak = findIndex(vfit, vpeak)[0] # index of vpeak on vfit
        iBfit = BL.linfit(vfit) # iBfit = linear fit of the baseline
        smallest_res = 1e9 # EMPIRICAL, typical value < 1e-3
        best_iWfit = None
        for deg in degrange:
            iWfit, res = WR.polynfit(vfit, deg) # iWfit = Degree n polynomial fit -- the working range
            _iIntersec, vIntersec = BL.intersection(vfit, iBfit, iWfit, _iVpeak) # index and value of the intersection of the two fits
            disjoint = len(vIntersec) == 0
            vbr = vIntersec[0] if not disjoint else np.nan
            #print(deg, res, _iIntersec, vIntersec, iBfit[_iIntersec])
            
            if res < smallest_res and not disjoint:
                smallest_res = res
                best_iWfit = iWfit
                bestvbr = vbr
                bestivbr = iBfit[_iIntersec][0]
                
        if "span" in pltItems:
            cdict = {BL: "c", WR: "y"} # color dict
            for segment in [BL, WR]:
                plt.axvspan(segment.xmin, segment.xmax, color=cdict[segment], alpha=0.05, zorder=2)
        if "linfit" in pltItems:
            _iBmin, _iBmax = findIndex(vfit, BL.xmin, vpeak+2) # indices of the endpoints of BL on vfit; Note: 2 volts right to vpeak; 2 volt ---> EMPIRICAL
            plt.plot(vfit[_iBmin:_iBmax], iBfit[_iBmin:_iBmax], c='c', ls = "-", lw=.8, zorder=4, label="") #f"{self.label} (linear)"
        if "polyfit" in pltItems:
            _iWmin, _iWmax = findIndex(vfit, WR.xmin, WR.xmax) # indices of the endpoints of WR on vfit
            plt.plot(vfit[_iWmin:_iWmax], best_iWfit[_iWmin:_iWmax], c='orange', ls = "-", lw=.8, zorder=4, label=f"{self.label} (Degree {deg})")
        if "intersect" in pltItems:
            plt.plot(bestvbr, bestivbr, 'ko', markersize=3)
        
        return bestvbr
    
    def calcVbr2(self):
        """
        Calculate the breakdown Voltage V_br by finding the leftmost peak on the curve of d/dV[ln(I)] vs. V. (Relative Derivative Method)

        Returns
        -------
        vbr : float
            Breakdown voltage, in volt.
            
        References
        ----------
        * (Relative derivative) \B. Bonanno et al. (2014). *IEEE Sens. J.* **14** (10) (2014) pp.3567–3578. https://doi.org/10.1109/JSEN.2014.2328623 
        * (Relative derivative) Hamamatsu MPPC® Technical Note (2022), https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/mppc_kapd9005e.pdf
        * (Summary of methods) \F. Nagy *et al.* (2017). *NIMPR, Section A* **849** pp.55-59. https://linkinghub.elsevier.com/retrieve/pii/S0168900217300025 
    
        """
        vpeak = self.calcPeak(order=1)
        vbr = vpeak - 0.18 
        # 0.18 V is the EMPIRICAL difference between the peak voltage found by the Relative Derivative Method and the breakdown voltage V_br.
        return vbr
    
    def calcVbr3(self):
        """
        Calculate the breakdown Voltage V_br by finding the leftmost peak on the curve of d^2/dV^2[ln(I)] vs. V. (Second Derivative Method)

        Returns
        -------
        vbr : float
            Breakdown voltage, in volt.
            
        References
        ----------
        * (Second derivative) \M. Simonetta *et al.* (2016). *NIMPR, Section A* **824** pp.146-147. https://doi.org/10.1016/j.nima.2015.11.023
        * (Summary of methods) \F. Nagy *et al.* (2017). *NIMPR, Section A* **849** pp.55-59. https://linkinghub.elsevier.com/retrieve/pii/S0168900217300025
        """
        vpeak = self.calcPeak(order=2)
        vbr = vpeak
        return vbr
  
class IVsegment():  
    def __init__(self, xax, yax, xmin, xmax):
        """
        Class constructor for the class ``IVsegment``, which segments a pair of arrays (xax, yax) into two parts, 
        fit one of them to a linear model and the other to a polynomial model, 
        and also allows to find the intersection of the two fits.
        
        This class possesses the attributes including all of its four parameters and two resulting sliced-off arrays.
        
        This class is specially designed for the method ``IVcurve.calcVbr1``.

        Parameters
        ----------
        xax : array-like, shape (n, )
            x-array to be segmented.
        yax : array-like, shape (n, )
            y-array to be segmented.
        xmin : float
            The smallest vlaue for the segmentation.
        xmax : float
            The largest vlaue for the segmentation.

        Returns
        -------
        None.

        """
        self.xmin, self.xmax = xmin, xmax #
        self.xall, self.yall = xax, yax # original arrays
        _imin, _imax = findIndex(xax, xmin, xmax) # the first and the last indices for selecting the IVsegment.
        self.xseg, self.yseg = xax[_imin:_imax], yax[_imin:_imax] # the resulting sliced-off arrays.
        
    def xfiner(self, num=1000):
        """
        Returns a 1-d array with the same endpoints but with a finer (or looser) spacing. 
        This can be used as an input of ``IVsegment.linfit`` or ``IVsegment.polynfit``.

        Parameters
        ----------
        num : int, optional
            Number of entries in the new array. The default is 1000.

        Returns
        -------
        samples
            Evenly spaced 1-d array with ``num`` entries.

        """
        samples = np.linspace(self.xall.min(), self.xall.max(), num)
        return samples
    
    def linfit(self, xfit):
        """
        Calculate the linear regression of the ``IVsegment``.

        Parameters
        ----------
        xfit : array-like, shape (num, )
            x-array for creating the linear-fit y-array.

        Returns
        -------
        yfit : array-like, shape (num, )
            The linear-fit for the y-array.
            
        References
        ----------
        .. [scipy.stats.linregress] https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
        """
        linreg = linregress(self.xseg, self.yseg)
        yfit = linreg[0] * xfit + linreg[1]
        return yfit
    
    def polynfit(self, xfit, deg):
        """
        Calculate the polynomial regression of the ``IVsegment``.

        Parameters
        ----------
        xfit : array-like, shape (num, )
            x-array for creating the polynomial fit.
        deg : int
            Degree of the fitting polynomial

        Returns
        -------
        yfit : array-like, shape (num, )
            The Degree ``deg`` polynomial fit of the y-array.
        residuals : float
            Residuals, the sum of squared residuals of the least squares fit.
            
        References
        ----------
        .. [numpy.polyfit] https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
       	pf = np.polyfit(self.xseg, self.yseg, deg, full=True)
        coef, residuals = pf[0], pf[1]
        yfit = np.polyval(coef, xfit)
        return yfit, residuals
        
    def intersection(self, x, y1, y2, index_low):
        """
        Find the intersection of two curves (each of which is an ``IVsegment``) defined by 2-d array.

        Parameters
        ----------
        x : array-like, shape (num, )
            The common x-array for both curve.
        y1 : array-like, shape (num, )
            The y-array for one curve.
        y2 : array-like, shape (num, )
            The y-array for the other curve.
        index_low : int 
            The lowest index to the find intersections.

        Returns
        -------
        index_intersec : int
            The index of the intersection on the shape-(num, ) array, e.g. ``x``.
        value_intersec : float
            The voltage value of the intersection.
            
        References
        ----------
        .. [StackOverflow] https://stackoverflow.com/a/28766902/22277122

        """
        ycompare = np.sign(y1[index_low:] - y2[index_low:])
        index_intersec = np.argwhere(np.diff(ycompare)).flatten()
        index_intersec += index_low
        value_intersec = x[index_intersec]
        return index_intersec, value_intersec

class Test():
    def __init__(self, Ts):
        self.Tli = Ts # Ts needs to be sorted in advance
        self.Tarray = np.array([float(T) for T in Ts])
        self.Vbrarray = np.empty((3,self.Tarray.shape[0]))
        self.RQarray = np.empty((2,self.Tarray.shape[0]))
    
    def folder(self):
        return "template_data"
    
    def filenames(self, T, region):
        return self.folder()  + '/T' + T + '_df_' + region + '.csv'
    
    def testIV(self, region, savePlot=False, printResult=False, placeLegend=True):
        """
        Plot the I-V curve from the template files and do the analysis.
    
        Parameters
        ----------
        region : {"r", "f"}
            Forward-biasing region (``"f"``) or reversed-biasing region (``"r"``).
        savePlot : bool, optional
            Whether to save the plot. The default is False.
        printResult : bool, optional
            Whether to print the result. The default is True.
        placeLegend : bool, optional
            Whether to place a legend. The default is True.
    
        Returns
        -------
        None.
    
        """
        if printResult==True:
            methods = "\n        | Method 1   Method 2   Method 3"
            headerdict = {"r": "Breakdown voltage V_br" + methods, "f": "Quenching resistance R_Q"}
            print("\nTemp.   | " + headerdict[region] + "\n" + "+"*40)
        
        yscaledict = {"r": "log", "f": "linear"}
        fnli = self.Tli
        IV = IVcurve(yscale=yscaledict[region]) 
        IV.setColors(fnli)
        
        for i, T in enumerate(fnli):
            label = str(int(T)) + " °C"
            fn = self.filenames(T, region) # generate the filename
            df = pd.read_csv(fn, sep=",") # open the file as a pandas.DataFrame
            
            # import the data
            if region == "r":
                iv = np.sort(df.to_numpy(), axis=0)
                x = iv[2:,0] # voltage
                y = iv[2:,1] # current
            if region == "f":
                x = df["Voltage"].to_numpy() # voltage
                y = df["Current"].to_numpy() # current
            
            IV.addPlot(x, y, color=IV.colorset[i], label=label)
    
            # analysis: find V_br for 
            if region == "r":
                IV.calcPeak(order=1, pltMode=0, pltItems=[]) # "vlines", "text"
                vbr1 = self.Vbrarray[0,i] = IV.calcVbr1(pltItems=[]) #"linfit", "polyfit", "intersect", "span"
                vbr2 = self.Vbrarray[1,i] = IV.calcVbr2()
                vbr3 = self.Vbrarray[2,i] = IV.calcVbr3()
                if printResult == True:
                    print(f"{int(fnli[i]):3d} °C  |  {vbr1:5.2f} V    {vbr2:5.2f} V    {vbr3:5.2f} V")
            if region == "f":
                Rq, Rqerr = IV.calcRQ(N=Nmicrocells, pltItems=[]) # "fit", "text"
                self.RQarray[:,i] = Rq, Rqerr
                if printResult == True:
                    print(f"{int(fnli[i]):3d} °C  |   ({Rq/1e3:7.2f} ± {Rqerr/1e3:.2f}) kΩ")
            if placeLegend == True:
                IV.addLegend(outright=True)
        
        if savePlot == True:
            IV.fig.savefig(TrialName + "-IV-" + region + ".png", dpi=200)

    def testVT(self, savePlot=False, printResult=True):
        """
        Plot the temperature dependence of breakdown voltage V_br and calculate the temperature coefficient ΔV/ΔT.
    
        Parameters
        ----------
        savePlot : bool, optional
            Whether to save the plot. The default is False.
        printResult : bool, optional
            Whether to print the result. The default is True.
    
        Returns
        -------
        None.
    
        """
        TvsVbr = np.vstack([self.Tarray, self.Vbrarray])
        fig, ax = plt.subplots(figsize=(7,4))
        prop = ax._get_lines.prop_cycler
        ax.set_xlabel("$T$ [°C]")
        ax.set_ylabel("$V_{br}$ [V]")
        xu = TvsVbr[0]
        xfi = np.linspace(xu.min(), xu.max(), 100)
        padli = []
        if printResult==True:
            print("")
        for i in [1,2,3]:
            color = next(prop)['color']
            yu = TvsVbr[i]
            p1, = ax.plot(xu, yu, "o", zorder=4, color=color)
            slope, interc, *other = linregress(xu, yu)
            yfi = xfi * slope + interc
            yufi = xu * slope + interc
            uncertainty = np.sqrt(1/(len(xu)-2) * ((yu-yufi)**2).sum() / ((xu - xu.mean())**2).sum())
            if printResult == True:
                print(f"Method {i}  : ΔV/ΔT = ({slope*1e3:.2f} ± {uncertainty*1e3:.2f}) mV/K")
            p2, = ax.plot(xfi, yfi, color=color)
            padli.append((p1,p2))
            
        zipdata = np.array([[TvsVbr[0][i], TvsVbr[1:].transpose()[i,j]] for i in range(TvsVbr.shape[1]) for j in range(3)])
        zipx =  zipdata[:,0]
        zipy =  zipdata[:,1]
        slope, interc, *other = linregress(zipx, zipy)   
        yfi = xfi * slope + interc
        zipyfi = zipx * slope + interc
        pc, = ax.plot(xfi, yfi, ls="--", lw=2.5, c="0.5", zorder=5)
        padli.append(pc)
        Vbrmethod = {1: "Intercept", 2: "Relative derivative", 3: "Second derivative"}
        labli = [Vbrmethod[i+1] for i in range(3)]
        labli.append("combined")
        uncertainty = np.sqrt(1/(len(zipx)-2) * ((zipy-zipyfi)**2).sum() / ((zipx - zipx.mean())**2).sum())
        if printResult==True:
            print(f"3 Methods : ΔV/ΔT = ({slope*1e3:.2f} ± {uncertainty*1e3:.2f}) mV/K")
        for j in range(TvsVbr.shape[1]):
            for i in [1,2,3]:
                yu = TvsVbr[i]
                #std = np.std(TvsVbr[1:,j], ddof=1)
                #ax.errorbar(xu, yu, yerr=std, fmt="", color='k', linestyle='', zorder=2)
        ax.legend(padli, labli, handler_map={tuple: HandlerTuple(ndivide=None)})
        if savePlot ==True:
            fig.savefig(TrialName + "-Vbr-T.png", dpi=200)
        return TvsVbr
       
    def testRT(self, savePlot=False, printResult=True):
        """
        Plot the temperature dependence of quenching resistance R_Q and calculate the temperature coefficient ΔR/ΔT.
    
        Parameters
        ----------
        savePlot : bool, optional
            Whether to save the plot. The default is False.
        printResult : bool, optional
            Whether to print the result. The default is True.
    
        Returns
        -------
        None.
    
        """
        TvsRQ = np.vstack([self.Tarray, self.RQarray])
        xr, yr, yerr = TvsRQ[0], TvsRQ[1], TvsRQ[2]
        xfi = np.linspace(xr.min(), xr.max(), 100)
        slope, interc, *other = linregress(xr, yr)   
        yfi = xfi * slope + interc
        yrfi = yr * slope + interc
        uncertainty = np.sqrt(1/(len(xr)-2) * ((yr-yrfi)**2).sum() / ((xr - xr.mean())**2).sum()) 
        fig, ax = plt.subplots(figsize=(7,4))
        ax.set_xlabel("$T$ [°C]")
        ax.set_ylabel("$R_Q$ [k$\Omega$]")
        ax.errorbar(xr, yr/1e3, yerr=yerr/1e3, marker="D", ls="", label="calc", zorder=2)
        ax.plot(xfi, yfi/1e3, label="fit",zorder=1)
        plt.legend()
        if printResult == True:
            print(f"\nΔR/ΔT =  ({slope:.2f} ± {uncertainty:.2f}) Ω/K")
        if savePlot ==True:
            fig.savefig(TrialName + "-RQ-T.png", dpi=200)
        return TvsRQ

class tTcurve(Curve):
	def __init__(self, xlabel="Time [s]", ylabel="Temperature [\u00b0C]", title="Temperature variation", yscale="linear"):
		Curve.__init__(self, xlabel, ylabel, title, yscale)

	def plot(self, data, cli=["g","b"]):
		self.ax.plot(data[0], data[1], color=cli[0], ls="-", lw=1, zorder=3, label="T6")
		self.ax.plot(data[0], data[2], color=cli[1], ls="-", lw=1, zorder=3, label="T7")
		
		
		
#==============================================================================

def findIndex(arr, *values):
    """
    Finds the index (indices) of the nearest value(s) in a numpy.array
    Idea from https://www.geeksforgeeks.org/find-the-nearest-value-and-the-index-of-numpy-array/
    Also from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    
    Parameters
    ----------
    arr : numpy.array
        An array.
    *values : tuple
        Values to be searched.

    Returns
    -------
    indices : list
        Index (indices) of the nearest value(s).

    """
    indices = []
    for val in values:
        diff_array = np.abs(arr - val)
        indices.append(diff_array.argmin())
    return indices

#==============================================================================
Nmicrocells = 14400 # number of microcells on our SiPM
sP = False # common switch of savePlot in the four test functions
pD = True # common switch of printResult in the four test functions
if sP == True:
    TrialName = input("Enter the common prefix of the filenames to save: ")
"""
# file numbers tha specifying the temperature
Ts         = ['00', '05', '10', '14', '20', '24', '30', '35', '40', '50', '60', '70']

mytest = Test(Ts)
mytest.testIV(region="f", savePlot=sP, printResult=pD)
mytest.testIV(region="r", savePlot=sP, printResult=pD)
mytest.testVT(savePlot=sP, printResult=pD)
mytest.testRT(savePlot=sP, printResult=pD)
"""
