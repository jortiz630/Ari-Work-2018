import numpy as np
import matplotlib as plt

tCold = 1.5e+04
rMin = 1
rMax = 30

def tempCyl(radius, tHot):
    return (((tCold - tHot)/(np.log(rMin)-np.log(rMax))) * (np.log(radius)-np.log(rMax))) + tHot
def tempSph(radius, tHot):
    return tHot - rMin*((tCold-tHot)/(rMax-rMin)) + ((1/radius)*((rMax*rMin*(tCold - tHot))/(rMax-rMin)))
#----------------------#
def kappaCyl(radius,tHot):
    return 5.6e-07*(tempCyl(radius,tHot))**(5/2)
def kappaSph(radius,tHot):
    return 5.6e-07*(tempSph(radius,tHot))**(5/2)
#----------------------#
def dTdr_sph(radius,tHot):
    return (radius**(-2))*((rMax*rMin*(tHot-tCold))/(rMax-rMin))*3.2408e-22
def dTdr_cyl(radius,tHot):
    return (radius**(-1))*((tCold-tHot)/(np.log(rMin)-np.log(rMax)))*3.2408e-22
#----------------------#
def heatCyl(radius,tHot):
    return kappaCyl(radius,tHot)*dTdr_cyl(radius,tHot)
def heatSph(radius,tHot):
    return kappaSph(radius,tHot)*dTdr_sph(radius,tHot)
#----------------------#
def L23f_temp(temperature, metallicity):
    L23f = np.zeros_like(temperature)
    tempB = 1.0e+06 + 1.5e+07*(metallicity)**(2.0/3.0)
    tempM = 1.5e+05
    tempR = 1.5e+04
    
    if (metallicity == 0.0):
        for i in np.arange(len(temperature)):
            alpha = -0.8
            tempB = 1.0e+06
            if ((temperature[i] >= tempR) and (temperature[i] <= tempB)):
                L23f[i] = 12.0*((temperature[i]/tempR)**(alpha))
                print i, temperature[i], (temperature[i]/tempR)
            elif (temperature[i] > tempB):
                L23f[i] = 12.0*((tempB/tempR)**(alpha))*((temperature[i]/tempB)**(1.0/3.0))
                print i, temperature[i]
            else:
                return 'error!!!'
        return L23f
    
    elif ((metallicity <= 1.0) and (metallicity > 0.0)):
        alpha = (1.0 + ((1.0/3.0)*np.log(metallicity)))
        for i in np.arange(len(temperature)):
            if ((temperature[i] >= tempR) and (temperature[i]<= tempM)):
                L23f[i] = 12.0*((temperature[i]/tempR)**(alpha))
                print i, temperature[i], (temperature[i]/tempR)
            elif ((temperature[i] > tempM) and (temperature[i] <= tempB)):
                L23f[i] = (12.0*(tempM/tempR)**(alpha))*(temperature[i]/tempM)**(-1.0)
                print i, temperature[i]
            elif (temperature[i] > tempB):
                L23f[i] = 12.0*((tempM/tempR)**(alpha))*((tempB/tempM)**(-1.0))*((temperature[i]/tempB)**(1.0/3.0))
            else:
                return 'error!!!'
            
        return L23f
    else:
        return 'too metal for me!!!'
#----------------------#
def read_coolingfunction(fname):
    Zg=[0.0,1.0E-03,1.0E-02,3.162278E-02,1.0E-01,0.31622777,1.0,3.1622777]
    Zstr=[]
    for Z in Zg:
        Zstr.append("{:.3f}".format(Z))
    logT,L0,L1,L2,L3,L4,L5,L6,L7=np.loadtxt(fname,unpack=True)
    plt.pyplot.plot(logT,L0,label=Zstr[0])
    plt.pyplot.plot(logT,L1,label=Zstr[1])
    plt.pyplot.plot(logT,L2,label=Zstr[2])
    plt.pyplot.plot(logT,L3,label=Zstr[3])
    plt.pyplot.plot(logT,L4,label=Zstr[4])
    plt.pyplot.plot(logT,L5,label=Zstr[5])
    plt.pyplot.plot(logT,L6,label=Zstr[6])
    plt.pyplot.plot(logT,L7,label=Zstr[7])
    plt.pyplot.xlabel("log(T)")
    plt.pyplot.ylabel(r"log($\Lambda$(T)) ergs/cm$^3$")
    plt.pyplot.legend()
    plt.pyplot.savefig("coolingfunc.pdf")
    plt.pyplot.show()
#----------------------#
def L23f_cyl(radius, tHot, metallicity):
    temperature = tempCyl(radius,tHot)
    L23f = np.zeros_like(temperature)
    tempB = 1.0e+06 + 1.5e+07*(metallicity)**(2.0/3.0)
    tempM = 1.5e+05
    tempR = 1.5e+04
    
    if (metallicity == 0.0):
        for i in np.arange(len(temperature)):
            alpha = -0.8
            tempB = 1.0e+06
            if ((temperature[i] >= tempR) and (temperature[i] <= tempB)):
                L23f[i] = 12.0*((temperature[i]/tempR)**(alpha))
                #print i, temperature[i], (temperature[i]/tempR)
            elif (temperature[i] > tempB):
                L23f[i] = 12.0*((tempB/tempR)**(alpha))*((temperature[i]/tempB)**(1.0/3.0))
                #print i, temperature[i]
            else:
                return 'error!!!'
        return L23f
    
    elif ((metallicity <= 1.0) and (metallicity > 0.0)):
        alpha = (1.0 + ((1.0/3.0)*np.log(metallicity)))
        for i in np.arange(len(temperature)):
            if ((temperature[i] >= tempR) and (temperature[i]<= tempM)):
                L23f[i] = 12.0*((temperature[i]/tempR)**(alpha))
                #print i, temperature[i], (temperature[i]/tempR)
            elif ((temperature[i] > tempM) and (temperature[i] <= tempB)):
                L23f[i] = (12.0*(tempM/tempR)**(alpha))*(temperature[i]/tempM)**(-1.0)
                #print i, temperature[i]
            elif (temperature[i] > tempB):
                L23f[i] = 12.0*((tempM/tempR)**(alpha))*((tempB/tempM)**(-1.0))*((temperature[i]/tempB)**(1.0/3.0))
                #print 12.0*((tempM/tempR)**(alpha))*((temperature[i]/tempB)**(1.0/3.0))
            else:
                return 'error!!!'
            
        return L23f
    else:
        return 'too metal for me!!!'
#----------------------#
def L23f_sph(radius, tHot, metallicity):
    
    temperature = tempSph(radius,tHot)
    L23f = np.zeros_like(temperature)
    tempB = 1.0e+06 + 1.5e+07*(metallicity)**(2.0/3.0)
    tempM = 1.5e+05
    tempR = 1.5e+04
    
    if (metallicity == 0.0):
        for i in np.arange(len(temperature)):
            alpha = -0.8
            tempB = 1.0e+06
            if ((temperature[i] >= tempR) and (temperature[i] <= tempB)):
                L23f[i] = 12.0*((temperature[i]/tempR)**(alpha))
                #print i, temperature[i], (temperature[i]/tempR)
            elif (temperature[i] > tempB):
                L23f[i] = 12.0*((tempB/tempR)**(alpha))*((temperature[i]/tempB)**(1.0/3.0))
                #print i, temperature[i]
            else:
                return 'error!!!'
        return L23f
    
    elif ((metallicity <= 1.0) and (metallicity > 0.0)):
        alpha = (1.0 + ((1.0/3.0)*np.log(metallicity)))
        for i in np.arange(len(temperature)):
            if ((temperature[i] >= tempR) and (temperature[i]<= tempM)):
                L23f[i] = 12.0*((temperature[i]/tempR)**(alpha))
                #print i, temperature[i], (temperature[i]/tempR)
            elif ((temperature[i] > tempM) and (temperature[i] <= tempB)):
                L23f[i] = (12.0*(tempM/tempR)**(alpha))*(temperature[i]/tempM)**(-1.0)
                #print i, temperature[i]
            elif (temperature[i] > tempB):
                L23f[i] = 12.0*((tempM/tempR)**(alpha))*((tempB/tempM)**(-1.0))*((temperature[i]/tempB)**(1.0/3.0))
                #print 12.0*((tempM/tempR)**(alpha))*((temperature[i]/tempB)**(1.0/3.0))
            else:
                return 'error!!!'
            
        return L23f
    else:
        return 'too metal for me!!!'
#----------------------#
def coolingFunctionCyl(radius, tHot, metallicity):
    return 1.0e-23 * L23f_cyl(radius, tHot, metallicity)
def coolingFunctionSph(radius, tHot, metallicity):
    return 1.0e-23 * L23f_sph(radius, tHot, metallicity)
#----------------------#
def radiativeHeatCyl(radius,tHot,metallicity):
    return ((10**(-3))**2) * coolingFunctionCyl(radius, tHot, metallicity)
def radiativeHeatSph(radius,tHot,metallicity):
    return ((10**(-3))**2) * coolingFunctionSph(radius, tHot, metallicity)
#----------------------#
def changingRadsCyl(radius,tHot,metallicity):
    return (((10**(3))/tempCyl(radius,tHot))**2) * coolingFunctionCyl(radius, tHot, metallicity)
def changingRadsSph(radius,tHot,metallicity):
    return (((10**(3))/tempSph(radius,tHot))**2) * coolingFunctionSph(radius, tHot, metallicity)
#----------------------#