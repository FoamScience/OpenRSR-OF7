""" Buckley-Lererett theory

Analytical solution to Buckley-Leverett equation with Brooks-Corey model
in 1D composite rock as it appears in the book (Multiphase fluid
flow in porous and fractured reservoirs) by YU-SHU WU (ISBN: 978-0-12-803848-2)

The main porpuse of this file is to provide theoretical values to validate
two-phase IMPES and coupled solvers of the OpenRSR project.

Author
------
name : Fadeli Mohamed Elwardi
email: elwardifadeli@gmail.com

"""

### Imports

import sys
import math
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt

### Classes

class rock():
    """rock struct storing rock-specific variables.
    
    Attributes
    ----------
    ID : int
        Rock id in global domain (used for plotting)
    length : float
        Rock length
    A : float
        Rock cross-section area
    porosity : float
        Effective porosity
    k : float, np.array
        Absolute permeability
    theta : float
        Angle between the 1D rock and the x-axis
    """

    def __init__(self, ID, length, A, porosity, k, theta):
        self.ID = ID
        self.length = length
        self.A = A
        self.porosity = porosity
        self.k = k
        self.theta = theta

class fluid():
    """fluid struct storing fluid properties.

    Attributes
    ----------
    name : str
        Fluid name
    mu : float
        Viscosity for Newtonian fluids
    rho : float
        Density for incompressible fluids
    """

    def __init__(self, name, mu, rho):
        self.name = name
        self.mu = mu
        self.rho = rho
        
class krProps():
    """ Kr operations for 2 fluid objs in a domain

    Attributes
    ----------
    fluid1 : fluid
        Wetting component
    fluid2 : fluid
        Non-Wetting component
    Smin1, n1, krMax1 : float
        Brooks Corey Coeffs for fluid1
    Smin2, n2, krMax2 : float
        Brooks Corey Coeffs for fluid2
    dS : Saturation step
    """

    def __init__(self, fluid1, fluid2, Smin1, Smin2, n1, n2, krMax1, krMax2, dS):
        self.fluid1 = fluid1
        self.fluid2 = fluid2
        self.Smin1 = Smin1
        self.n1 = n1
        self.krMax1 = krMax1
        self.Smin2 = Smin2
        self.n2 = n2
        self.krMax2 = krMax2
        self.dS = dS

    def kr1(self, S):
        """Kr curves for wetting phase"""
        return self.krMax1*(( S+eps-self.Smin1)/(1-self.Smin1-self.Smin2))**self.n1

    def kr2(self, S):
        """Kr curves for nonwetting phase"""
        return self.krMax2*((eps+1-S-self.Smin2)/(1-self.Smin1-self.Smin2))**self.n2
        
    def mob(self, S):
        """Mobility ratio"""
        return (self.kr2(S)/self.fluid2.mu)*(self.fluid1.mu/self.kr1(S))

### Helper Functions and global variables

#### Machine epsilon
eps = sys.float_info.epsilon

#### Time step
gdt = 86400*0.1

#### Starting and ending time
t0 = gdt
tf = 10*gdt

#### Gravitational acceleration
g = 9.81

#### Injected flowrate in domain 1
qt = 1e-5

#### Reference pressure (chosen at producer)
Pref = 0

def fw(S, krObj, rk, qt):
    """Calculate fractional flow array in a domain with injected flowrate qt"""

    dRho = krObj.fluid1.rho - krObj.fluid2.rho
    return 1/(1+krObj.mob(S)) \
            - rk.A * rk.k * krObj.kr1(S) / ( qt * krObj.fluid1.mu) * \
            dRho * g * math.sin(rk.theta)/(1+krObj.mob(S))

def dfw(S, krObj, rk, qt):
    """Calculate fractional flow derivative array in a domain"""

    dRho = krObj.fluid1.rho - krObj.fluid2.rho
    return (( krObj.n1 * krObj.mob(S) / (S-krObj.Smin1)) \
             +(krObj.n2 * krObj.mob(S)) / (1-S-krObj.Smin2)) \
           / ( 1 + krObj.mob(S) )**2 \
           + rk.A * rk.k * krObj.kr2(S) * krObj.n2 * dRho * g \
             * math.sin(rk.theta) \
           / ( (1 - S - krObj.Smin2) * qt * krObj.fluid2.mu \
             * (1 + krObj.mob(S) ) ) \
           + rk.A * rk.k * krObj.kr2(S) * dRho * g * math.sin(rk.theta) \
           / ( qt * krObj.fluid2.mu * (1 + krObj.mob(S))**2) \
           * (((krObj.n1 * krObj.mob(S) / (S - krObj.Smin1)) \
           + (krObj.n2 * krObj.mob(S)) /(1 - S - krObj.Smin2)) \
           / (1 + krObj.mob(S))**2)

def frontSaturation(S, fw, dfw):
    """Estimate front saturation properties using Welge method

    Returns
    -------
    index : int
        index of advance front saturation in saturation array
    Swf : float
        Advancing front saturation
    dfwf : float
        dfw at front saturation
    bf : float
        Intercept value for Weldge's tangent line
    """

    # Initiate returns:
    index, Swf, dfwf, bf = np.nan, np.nan, np.nan, np.nan
    ns = len(S)
    df = 1/(S[ns-1] - S[0])

    for i in range(1, ns-1):
        dfds = fw[i] / (S[i] - S[0])
        if ((dfds < dfw[i-1]) and (dfds > dfw[i+1]) and (dfds >= df)):
            index = i
            Swf = S[i]
            dfwf = dfds
            bf = fw[i] - dfwf * S[i]
            break

    if np.isnan(index):
        raise Exception("dS in krProps is too big!!")

    return (index, Swf, dfwf, bf)

def satBeforeInterfaceReach(S, fwArray, dfwArray, krObj, rk, qt, t0):
    """Plot and write saturation profile when time < tf

    Returns
    -------
    tff : float
        Time to reach interface of domains
    index: int
        front saturation index in saturation array
    x : np.array
        Positions through domains to calculate saturations at
    Sw : np.array
        Saturation values at picked position and time values
    """

    ### Frontal water saturation position
    (index, Swf, dfwf, bf) = frontSaturation(S, fwArray, dfwArray)

    ### Water saturation profile

    krw_swf = krObj.kr1(Swf)
    kro_swf = krObj.kr2(Swf)
    fw_swf = fw(Swf, krObj, rk, qt)
    dfw_swf = dfw(Swf, krObj, rk, qt)

    tff = rk.A * rk.porosity * rk.length / ( qt * dfw_swf)

    t = np.arange(t0, tff if tff <= tf else tf, step=gdt)
    ndt = len(t)

    St = np.flip(S)
    dfwt = np.flip(dfwArray)

    Xsw = np.zeros(shape=(len(S)-index+1, ndt))
    XX = [[] for i in range(len(t))]
    SS = [[] for i in range(len(t))]
    for j in range(ndt):
        for k in range(len(S)-index+1):
            Xsw[k,j] = qt * t[j] / ( rk.A * rk.porosity ) * dfwt[k]
    plt.figure(1)
    plt.ylabel('$S_w$')
    plt.xlabel('x [m]')
    plt.xlim([0, rk.length*1.01])
    plt.ylim([0, 1])
    for j in range(len(t)):
        x = np.append([0, Xsw[0,j]], Xsw[:,j])
        x = np.append(x, [Xsw[-1,j], Xsw[-1,j]])
        x = np.append(x, [Xsw[-1,j], rk.length])
        Sw = np.append([1-krObj.Smin2, 1-krObj.Smin2], St[:len(S)-index+1])
        Sw = np.append(Sw, [Swf, krObj.Smin1])
        Sw = np.append(Sw, [krObj.Smin1, krObj.Smin1])
        plt.plot(x, Sw, label="{:.2} days".format(t[j]/86400.0))
        df = pd.DataFrame({'x': x, 'S':Sw})
        df.to_csv("water.alpha/water.alpha-{:.2f}-days.csv".format(t[j]/86400.0), sep=" ",
                header=None, index=None)
        XX[j] = x
        SS[j] = Sw
    plt.legend()
    plt.title("Saturation profile in domain {}".format(rk.ID))

    return (tff, index, XX, SS)

def satAfterInterfaceReach(S, fwArray, dfwArray, krObj, rk, qt, t0, index):
    """Calculate saturation profile when time >= tf

    Returns
    -------
    return (SwAvg, fwf, Xsw)
    SwAvg : float
        Average saturation in domain 1 (For oil recovery calculation)
    fwf: float
        Fractional flow at font saturation
    Xsw : np.array
        Positions where chosen saturation values are
    """
    t = np.flip(np.arange(tf, ti, step=-gdt))
    t = np.append([ti], t)
    ndt=len(t)
    dfwf = rk.A * rk.porosity * rk.length / qt / t
    interpFunc = interp(dfwArray[index:], S[index:])
    Swf = interpFunc(dfwf)
    #krwf = krObj.kr1(Swf)
    #krof = krObj.kr2(Swf)
    fwf = fw(Swf, krObj, rk, qt)
    SwAvg = Swf + (1-fwf) / dfwf
    
    dfwt = np.flip(dfwArray)

    Xsw = np.zeros(shape=(len(S)-index+1, ndt))
    for j in range(ndt):
        for k in range(len(S)-index+1):
            Xsw[k,j] = qt * t[j] / ( rk.A * rk.porosity ) * dfwt[k]
    return (SwAvg, fwf, Xsw)


def pressureProfile(Pref, Xsw, S, rk1, krObj1, rk2, krObj2, time):
    """Calculate Pressure profile in a domains

    """


    dx = (rk1.length + rk2.length)/100
    x = np.arange(rk1.length+rk2.length, 0, -dx)

    P = np.linspace(Pref, Pref, len(x))

    plt.figure(4)
    interpFunc = interp(Xsw, S)
    for j in range(1, len(x)):
        if (x[j] < rk1.length):
            swp = interpFunc(x[j])
            fwp = fw(swp, krObj1, rk1, qt)
            Gp = krObj1.fluid2.rho * g * math.sin(rk1.theta)
            Tp = rk1.A * rk1.k * krObj1.kr2(swp) / krObj1.fluid2.mu
            P[j] = dx * (qt * (1 - fwp) / Tp - Gp) + P[j-1]
        elif (x[j] >= rk1.length and max(Xsw) > rk1.length):
            swp = interpFunc(x[j])
            fwp = fw(swp, krObj2, rk2, qt)
            Gp = krObj2.fluid2.rho * g * math.sin(rk2.theta)
            Tp = rk2.A * rk2.k * krObj2.kr2(swp) / krObj2.fluid2.mu
            P[j] = dx * (qt * (1 - fwp) / Tp - Gp) + P[j-1]
        else:
            P[j] = Pref

    plt.plot(x, P)

    df = pd.DataFrame({'x': x, 'P':P})
    df.to_csv("p/p-{:.2f}-days.csv".format(time/86400.0), sep=" ",
            header=None, index=None)

    return (x, P)



if __name__ == "__main__":

    ### Fluid properties
    water = fluid("Water", 1e-3, 1e3)
    oil = fluid("Oil", 5e-3, 0.8e-3)

    ### Paramerter initialization for domain 1
    rk1 = rock(1, 6.0, 1.0, 0.3, 1e-14, 0.0)
    krO1 = krProps(water, oil, 0.2, 0.2, 1.5, 2.5, 0.8, 0.8, 1e-3)
    S1 = np.arange(krO1.Smin1+eps, 1-krO1.Smin1-eps, step=krO1.dS)

    ### Paramerter initialization for domain 2
    rk2 = rock(2, 6.0, 1.0, 0.3, 1e-14, 0.0)
    krO2 = krProps(water, oil, 0.2, 0.2, 2.5, 1.5, 0.75, 0.75, 1e-3)
    S2 = np.arange(krO2.Smin1+eps, 1-krO2.Smin1-eps, step=krO1.dS)

    ### Fractional flow functions

    krw1 = krO1.kr1(S1)
    kro1 = krO1.kr2(S1)
    fw1 = fw(S1, krO1, rk1, qt)
    dfw1 = dfw(S1, krO1, rk1, qt)

    krw2 = krO2.kr1(S1)
    kro2 = krO2.kr2(S1)
    fw2 = fw(S2, krO2, rk2, qt)
    dfw2 = dfw(S2, krO2, rk2, qt)

    plt.figure(2)
    plt.ylabel('$k_r$')
    plt.xlabel('$S_w$')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.plot(S1, krw1, label="{} Kr in domain {}".format(water.name, rk1.ID))
    plt.plot(S1, kro1, label="{} Kr in domain {}".format(oil.name, rk1.ID))
    plt.plot(S2, krw2, "--", label="{} Kr in domain {}".format(water.name, rk2.ID))
    plt.plot(S2, kro2, "--", label="{} Kr in domain {}".format(oil.name, rk2.ID))
    plt.legend()
    plt.title("Brooks-Corey relative permeability curves")
    #plt.show()

    plt.figure(3)
    plt.ylabel('$f_w$ and $dfw$')
    plt.xlabel('$S_w$')
    plt.xlim([0, 1])
    plt.plot(S1, fw1, label="fw in domain {}".format(rk1.ID))
    plt.plot(S1, dfw1,"--", label="dfw in domain {}".format(rk1.ID))
    plt.plot(S2, fw2, label="fw in domain {}".format(rk2.ID))
    plt.plot(S2, dfw2, "--", label="dfw in domain {}".format(rk2.ID))
    plt.legend()
    plt.title("Fractional flow curves")
    #plt.show()

    ### Saturation profiles in domain 1

    #### Before raching domains interface
    (ti, index, xr1, sr1) = satBeforeInterfaceReach(S1, fw1, dfw1, krO1,
            rk1, qt, t0)
    for i in range(len(np.arange(t0, ti, step=gdt))):
        xv = np.append(xr1[i], [rk1.length+rk2.length])
        Sv = np.append(sr1[i], [krO2.Smin1])
        pressureProfile(Pref, xv, Sv, rk1, krO1, rk2, krO2, t0 + gdt*i)

    if (ti > tf):
        print("Flow didn't reach domain 2")
        plt.show()
        sys.exit()

    #### After reaching domain2
    (Sw1Avg, fwf1, xr2) = satAfterInterfaceReach(S1, fw1, dfw1, krO1, rk1,
            qt, ti, index)

    #for i in range(len(np.arange(t0, ti, step=gdt))):
    #    xs = xr2[:, i]
    #    xv = xs[xs < rk1.length]
    #    pressureProfile(Pref, np.flip(xv), S1[:len(xv)], rk1, krO1, rk2,
    #            krO2, t0 + gdt*i)

    ### Saturation profiles in domain 2
    #t2 = np.linspace(ti, tf, ndt2)
    t2 = np.flip(np.arange(tf, ti, step=-0.1*86400))
    t2 = np.append([ti], t2)
    ndt2 = len(t2)
    interFunc = interp(fw2, S2, kind='cubic')
    Sw2f = interFunc(fwf1)
    #Sw2f[2] = 0.511
    Qt = qt * t2
    V1 = rk1.length * rk1.A * rk1.porosity * (Sw1Avg - krO1.Smin1)
    V2 = Qt - V1
    if (abs(V2[0]/Qt[0]) <= 1e-3):
        V2[0] = 0

    Sw2_k = np.ndarray(shape=(ndt2,1), dtype=np.ndarray)
    Xsw2_k = np.ndarray(shape=(ndt2, 1), dtype=np.ndarray)
    Sw2_k[0,0] = np.array([Sw2f[0], krO2.Smin1])
    Xsw2_k[0,0] = np.array([rk1.length, rk1.length])

    index2 = [0 for i in range(ndt2)]
    Sw2fk = np.ndarray(shape=(ndt2,1), dtype=np.ndarray)
    Xswf2k = np.ndarray(shape=(ndt2,1), dtype=np.ndarray)
    for i in range(1, ndt2):
        Sw2_k1 = np.arange(Sw2f[i]-krO2.dS, Sw2f[0], step=-krO2.dS)
        Sw2_k2 = np.arange(Sw2f[0]-krO2.dS, krO2.Smin1, step=-krO2.dS)
        Sw2_k[i,0] = np.concatenate([[Sw2f[i]], Sw2_k1, Sw2_k2]).ravel()

        krw2_k1 = krO2.kr1(Sw2_k1)
        kro2_k1 = krO2.kr2(Sw2_k1)
        fw2_k1 = fw(Sw2_k1, krO2, rk2, qt)
        dfw2_k1 = dfw(Sw2_k1, krO2, rk2, qt)
        fw1_k1 = fw2_k1
        interpFunc = interp(fw1, S1)
        Sw1_k1 = interpFunc(fw1_k1)
        krw1_k1 = krO1.kr1(Sw1_k1)
        kro1_k1 = krO1.kr2(Sw1_k1)
        fw1_k1 = fw(Sw1_k1, krO1, rk1, qt)
        dfw1_k1 = dfw(Sw1_k1, krO1, rk1, qt)
        ts = rk1.A * rk1.porosity * rk1.length / qt / dfw1_k1
        Xsw2_k1 = rk1.length + qt / rk2.A / rk2.porosity * dfw2_k1 * (t2[i]-ts)
        krw2_k2 = krO2.kr1(Sw2_k2)
        kro2_k2 = krO2.kr2(Sw2_k2)
        fw2_k2 = fw(Sw2_k2, krO2, rk2, qt)
        dfw2_k2 = dfw(Sw2_k2, krO2, rk2, qt)
        Xsw2_k2 = rk1.length + qt / rk2.A / rk2.porosity * dfw2_k2*(t2[i]-ti)
        Xsw2_k[i,0] = np.concatenate([[rk1.length], Xsw2_k1, Xsw2_k2]).ravel()
        dXsw2 = np.array(Xsw2_k[i,0][1:-1])-Xsw2_k[i,0][0:-2]
        dSwk2 = np.array(Sw2_k[i,0][1:-1])-krO2.Smin1

        for j in range(len(S2)):
            W2 = rk2.A * rk2.porosity * sum(dXsw2[0:j+1]*dSwk2[0:j+1])
            if (W2 >= V2[i]):
                index2[i] = j
                Sw2fk[i] = Sw2_k[i,0][j-1]
                Xswf2k[i] = Xsw2_k[i,0][j-1]
                break

        Xsw2_k[i,0] = Xsw2_k[i,0][0:index2[i]]
        Sw2_k[i,0] = Sw2_k[i,0][0:index2[i]]

    plt.figure(1)
    L = rk1.length + rk2.length
    plt.ylabel('$S_w$')
    plt.xlabel('x [m]')
    plt.xlim([0, L*1.01])
    plt.plot([rk1.length, rk1.length], [0,1], "-.c", alpha=0.5)
    plt.text(rk1.length*1.05, 0.95, s="Domain 2")
    plt.text(rk1.length*0.95, 0.95, s="Domain 1", ha='right')
    for j in range(ndt2):
        Xfrom1 = [xr2[i,j] for i in range(len(xr2[:,j])) if xr2[i,j] <= rk1.length]
        Xsw2_k[j,0] = np.append(Xfrom1, Xsw2_k[j,0])
        Xsw2_k[j,0] = np.append(Xsw2_k[j,0], [Xsw2_k[j,0][-1], Xsw2_k[j,0][-1]])
        Xsw2_k[j,0] = np.append(Xsw2_k[j,0], [Xsw2_k[j,0][-1], L])
        Sw2_k[j,0] = np.append(np.flip(S1)[:len(Xfrom1)], Sw2_k[j,0])
        Sw2_k[j,0] = np.append(Sw2_k[j,0], [Sw2_k[j,0][-1], krO2.Smin1])
        Sw2_k[j,0] = np.append(Sw2_k[j,0], [krO2.Smin1, krO2.Smin1])
        plt.plot(Xsw2_k[j,0], Sw2_k[j,0], label="{:.2} days".format(t2[j]/86400.0))
        df = pd.DataFrame({'x': Xsw2_k[j,0], 'S':Sw2_k[j,0]})
        df.to_csv("water.alpha/water.alpha-{:.2f}-days.csv".format(t2[j]/86400.0),
                sep=" ", header=None, index=None)

    plt.plot([0, L], [Sw2f[0], Sw2f[0]], "--k", alpha=0.7)
    plt.text(L*0.98, Sw2f[0]*1.02, s='$S_{wf}$')
    plt.legend()
    plt.title("Saturation profile through the domains")


    #### Pressure profile while front is in domain 2
    for i in range(ndt2):
        pressureProfile(Pref, Xsw2_k[i,0], Sw2_k[i, 0], rk1, krO1, rk2, krO2, t2[i])
    plt.figure(4)
    plt.ylabel('P [Pa]')
    plt.xlabel('x [m]')
    plt.xlim([0, L*1.01])
    plt.title("Pressure profile through domains")
    plt.axvline(x=rk1.length, ymin=0, ymax=1, c="c", ls='-.', alpha=0.5)
    plt.text(rk1.length*1.05, 0.95, s="Domain 2")
    plt.text(rk1.length*0.95, 0.95, s="Domain 1", ha='right')
    plt.show()
