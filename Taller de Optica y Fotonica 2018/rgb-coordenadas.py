#import pylab
import matplotlib.pylab as plt
#import matplotlib.pyplot as plt # tambien funciona
import numpy as np 
import colour
from random import random, seed
seed(20180501)
from iminuit import Minuit, describe, Struct
from iminuit.util import make_func_code # lo utiliza Chi2Functor
from math import pi, exp, sqrt
#import colour.colorimetry as colorimetry

"""
Conviete valores RGB cualquiera, dado como constante, a xy y lo grafica como un punto dentro chromaticity diagram standar 
"""
from colour.plotting import (
    colour_plotting_defaults, chromaticity_diagram_plot_CIE1931,render,single_spd_plot)

# Defining a sample spectral power distribution data.
sample_spd_data = {380+i*5: random() for i in range(80)}
spd = colour.SpectralPowerDistribution(sample_spd_data)
single_spd_plot(spd)
cmfs = colour.STANDARD_OBSERVERS_CMFS['CIE 1931 2 Degree Standard Observer']
XYZ = colour.spectral_to_XYZ(spd, cmfs)
#ch1=colour.XYZ_to_spectral(XYZ,Method='Meng 2015')
ch1=colour.XYZ_to_spectral(XYZ,Method='Smits 1999')
single_spd_plot(ch1)
#spd_rec=ch1.values
y=ch1.values
x=[360+i*5 for i in range(95)]
#for i, v in enumerate(spd_rec):
#    x=360+i*5
#    y=v
#this is very useful if you want to build a generic cost functor
#this is actually how probfit is implemented
#x=sorted(sample_spd_data.keys())
#y=sorted(sample_spd_data.values())
class Chi2Functor:
    def __init__(self,f,x,y):
        self.f = f
        self.x = x
        self.y = y
        f_sig = describe(f)
        #this is how you fake function 
        #signature dynamically
        self.func_code = make_func_code(f_sig[1:])#docking off independent variable
        self.func_defaults = None #this keeps np.vectorize happy
    def __call__(self,*arg):
        #notice that it accept variable length
        #positional arguments
        chi2 = sum((y-self.f(x,*arg))**2 for x,y in zip(self.x,self.y))
        return chi2

def gaussian3(x,p1,p2,p3):
    g=(p1*exp(-(x-440)**2)+p2*exp(-(x-530)**2)+p3*exp(-(x-620)**2))
    return g/sqrt(pi)

comp3=Chi2Functor(gaussian3,x,y)
describe(comp3)
m=Minuit(comp3,p1=1,p2=1,p3=1)
m.migrad();
print(m.values)
sample_spd_fit = {i:gaussian3(i,m.values['p1'],m.values['p2'],m.values['p3']) for i in x }

spd_fit = colour.SpectralPowerDistribution(sample_spd_fit)
XYZf = colour.spectral_to_XYZ(spd_fit, cmfs)
xy = colour.XYZ_to_xy(XYZ) # Conversion to correlated colour temperature in K. 
xyf =  colour.XYZ_to_xy(XYZf)
print("XYZ: ",XYZ)
print(xy)
print("XYZf: ",XYZf)
print(xyf)
#print(colour.delta_E(XYZ,XYZf))

# Plotting the *CIE 1931 Chromaticity Diagram*.
# The argument *standalone=False* is passed so that the plot doesn't get displayed
# and can be used as a basis for other plots.
chromaticity_diagram_plot_CIE1931(standalone=False)

# Plotting the *xy* chromaticity coordinates.
x, y = xy
xf,yf=xyf
print(sqrt((x-xf)**2+(y-yf)**2))

#pylab.plot(x, y, 'o-', color='black')
plt.plot(x, y, 'o-', color='black')
plt.plot(xf, yf, 'o-', color='blue')
plt.show()

#illuminant = colour.ILLUMINANTS_RELATIVE_SPDS['D65']
# Calculating the sample spectral power distribution *CIE XYZ* tristimulus values.
#
#XYZ = colour.spectral_to_XYZ(spd, cmfs, illuminant)
#spd_rec=colour.XYZ_to_spectral(XYZ)
#colour.REFLECTANCE_RECOVERY_METHODS.keys()

