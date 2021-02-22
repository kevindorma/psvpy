#!/usr/bin/env python
# coding: utf-8

# # Calculation Template
# ## Client: INTERNAL
# ---
# ## Project: PSV tool example
# ## Calc: 2020-CALC-PSV-001
# ## By: K. Dorma
# ## Date: December, 2020
# ---
# ## Authentication
# > Stamp, Permit
# ---
# ## Revision History
# |Revision | Date | Description | By | Reviewer|
# | :-------| :----|:------------|:---|:--------|
# |    1.0  | Dec. 2020 | Demo code | KCD |  |
# |    2.0  | feb 13 2020   | Python    | KCD     |    |
# 
# ---

# In[5]:


import numpy as np
import pandas as pd
import math
import scipy as sp
from scipy import interpolate

def waterTsat(P_kPa):
    # equation for saturated steam, 1/T = A + B*ln P + C/ln P
    # I fit this myself, 100 to 20000 kPaa
    # pressure is in kPaa
    # we need to replace this with IAPWS IFC97 formulas

    lnP = math.log(P_kPa)
 
    # constants for fit of 1/T = A + B*ln P + C/ln P
    tA = 0.00379302
    tB = -0.000220828
    tC = -0.000425693

    invT = tA + tB*lnP + tC/lnP
    return((1.0/invT) - 273.15)


def waterPsat(T_c):
    # equation for saturated steam, ln P_Pa = A + B/Tk + C*Tk + D*ln(Tk)
    # I fit this myself, 100 to 20000 kPaa
    # pressure is in kPaa

    T_k = T_c + 273.15

    pA = 116.6408494
    pB = -8572.035364
    pC = 0.013736471
    pD = -14.73621925
    lnPsat = pA + pB/T_k + pC*T_k + pD*math.log(T_k)
    return(math.exp(lnPsat)/1000.0)


def getKsh(PkPa, State):
    # State is either a temperature in degC or a string
    # if this is a temperature, then look up the value in the table
    # if not, just return 1.0 because this is saturated
    # throw an error if the tempeature is below saturation
    # tables for superheat factors. I had a larger table but it threw weird errors. Smaller seems to be better.
    
    Ksh_table = np.full((26,10), 0.0)    # pressure in rows, temperature in columns

    # julia code for Ksh temperature values in degC
    Ksh_tC = np.array([93.33333333,
		148.8888889,
		204.4444444,
		260.0,
		315.5555556,
		371.1111111,
		426.6666667,
		482.2222222,
		537.7777778,
		565.5555556]);
    
        # julia code for Ksh pressure values in kPa   
    Ksh_pkPa = np.array([137.8951817	,
		344.7379543	,
		689.4759087	,
		1034.213863	,
		1378.951817	,
		1723.689772	,
		2068.427726	,
		2413.16568	,
		2757.903635	,
		3102.641589	,
		3447.379543	,
		3792.117498	,
		4136.855452	,
		4826.331361	,
		5515.807269	,
		6205.283178	,
		6894.759087	,
		7584.234995	,
		8273.710904	,
		8963.186813	,
		9652.662721	,
		10342.13863	,
		12065.8284	,
		13789.51817	,
		17236.89772	,
		20684.27726	]);
    
    
    # Julia code for Ksh, rows are pressure, columns are temperature    
    # subcooled values are denoted with KSH = 1.
    # Interpolation near the saturation temperature could give artificially low value for Ksh
    # a better method might be to replace the value 1 for the nearest subcooled temperature
    # with the value that gives 1 when interpolated to the saturation tempeature.
    # Clumsy, but it should work
    Ksh_table = np.array([[	1	,	0.99455814	,	0.987	,	0.93	,	0.882	,	0.841	,	0.805	,	0.774	,	0.745	,	0.732	],
	[	1	,	0.997925224	,	0.987	,	0.93	,	0.882	,	0.841	,	0.805	,	0.774	,	0.745	,	0.732	],
	[	1	,	1	,	0.998	,	0.935	,	0.885	,	0.843	,	0.807	,	0.775	,	0.746	,	0.733	],
	[	1	,	1	,	0.984	,	0.94	,	0.888	,	0.846	,	0.808	,	0.776	,	0.747	,	0.733	],
	[	1	,	1	,	0.979	,	0.945	,	0.892	,	0.848	,	0.81	,	0.777	,	0.748	,	0.734	],
	[	1	,	1	,	1	,	0.951	,	0.895	,	0.85	,	0.812	,	0.778	,	0.749	,	0.735	],
	[	1	,	1	,	1	,	0.957	,	0.898	,	0.852	,	0.813	,	0.78	,	0.75	,	0.736	],
	[	1	,	1	,	1	,	0.963	,	0.902	,	0.854	,	0.815	,	0.781	,	0.75	,	0.736	],
	[	1	,	1	,	1	,	0.963	,	0.906	,	0.857	,	0.816	,	0.782	,	0.751	,	0.737	],
	[	1	,	1	,	1	,	0.961	,	0.909	,	0.859	,	0.818	,	0.783	,	0.752	,	0.738	],
	[	1	,	1	,	1	,	0.961	,	0.914	,	0.862	,	0.82	,	0.784	,	0.753	,	0.739	],
	[	1	,	1	,	1	,	0.962	,	0.918	,	0.864	,	0.822	,	0.785	,	0.754	,	0.74	],
	[	1	,	1	,	1	,	0.964	,	0.922	,	0.867	,	0.823	,	0.787	,	0.755	,	0.74	],
	[	1	,	1	,	1	,	1	,	0.931	,	0.872	,	0.827	,	0.789	,	0.757	,	0.742	],
	[	1	,	1	,	1	,	1	,	0.942	,	0.878	,	0.83	,	0.792	,	0.759	,	0.744	],
	[	1	,	1	,	1	,	1	,	0.953	,	0.883	,	0.834	,	0.794	,	0.76	,	0.745	],
	[	1	,	1	,	1	,	1	,	0.959	,	0.89	,	0.838	,	0.797	,	0.762	,	0.747	],
	[	1	,	1	,	1	,	1	,	0.962	,	0.896	,	0.842	,	0.8	,	0.764	,	0.749	],
	[	1	,	1	,	1	,	1	,	0.966	,	0.903	,	0.846	,	0.802	,	0.766	,	0.75	],
	[	1	,	1	,	1	,	1	,	0.973	,	0.91	,	0.85	,	0.805	,	0.768	,	0.752	],
	[	1	,	1	,	1	,	1	,	0.982	,	0.918	,	0.854	,	0.808	,	0.77	,	0.754	],
	[	1	,	1	,	1	,	1	,	0.993	,	0.926	,	0.859	,	0.811	,	0.772	,	0.755	],
	[	1	,	1	,	1	,	1	,	1	,	0.94	,	0.862	,	0.81	,	0.77	,	0.752	],
	[	1	,	1	,	1	,	1	,	1	,	0.952	,	0.861	,	0.805	,	0.762	,	0.744	],
	[	1	,	1	,	1	,	1	,	1	,	0.951	,	0.852	,	0.787	,	0.74	,	0.721	],
	[	1	,	1	,	1	,	1	,	1	,	1	,	0.831	,	0.753	,	0.704	,	0.684	]])

# we have tables
        # we will use linear interpolation with P and T. 
        # Using ln P and 1/T might be more robust and needs to be investigated.

    # I need to find an interpolation routine
#    linear_Ksh = reshape(Ksh_table,(10*26,1));
#    Ksh_grid = GridInterpolations.RectangleGrid(Ksh_pkPa, Ksh_tC);  	# rectangular grid
#>>> xx, yy = np.meshgrid(x, y)
#>>> z = np.sin(xx**2+yy**2)
# this is the interpolation function
    returnFnct = sp.interpolate.interp2d(Ksh_tC, Ksh_pkPa, Ksh_table, kind='linear')

        
    Ksh = 1.0                # default value
    if (isinstance(State,float)):
        returnVal = returnFnct(State, PkPa)
        Ksh = returnVal[0]
        # check if we are subcooled
        pSat = waterPsat(State)
        if (pSat < PkPa):
            raise Exception("Temperature is in subcooled region")

    return (Ksh)
# that was easy


def PSVareaOrifice(letter):
# from the PSV designation letter, output the PSV area in MM2
# create the table
    data = {"Designation": ["D","E","F","G","H","J","K","L","M","N","P","Q","R","T"],
        "typicalFlanges": ["1.5D2", "1.5E2", "1.5F3", "1.5G3", "2H3", "3J4", "3K4", "4L6", "4M6", "4N6", "4P6", "6Q8", "6R10", "8T10"],
        "areaIN2": [0.11, 0.20, 0.31, 0.50, 0.79, 1.29, 1.84, 2.85, 3.60, 4.34, 6.38, 11.05, 16.00, 26.00],
        "areaMM2": [70.9676, 126.4514, 198.0641, 324.5155, 506.4506, 830.3209, 1185.804, 1840.641, 2322.576, 2799.994, 4116.121, 7129.018, 10322.56, 16774.16] }

    psvTable = pd.DataFrame(data, columns = ["Designation", "typicalFlanges", "areaIN2", "areaMM2"])


    returnVal = psvTable[psvTable['Designation'].str.contains(letter)]
    return(returnVal.iloc[0]['areaMM2'])





def PSVdesignationOrifice(myAreaMM2):
# from the aream, output the designation for the next larger orifice
# create the table
    data = {"Designation": ["D","E","F","G","H","J","K","L","M","N","P","Q","R","T"],
        "typicalFlanges": ["1.5D2", "1.5E2", "1.5F3", "1.5G3", "2H3", "3J4", "3K4", "4L6", "4M6", "4N6", "4P6", "6Q8", "6R10", "8T10"],
        "areaIN2": [0.11, 0.20, 0.31, 0.50, 0.79, 1.29, 1.84, 2.85, 3.60, 4.34, 6.38, 11.05, 16.00, 26.00],
        "areaMM2": [70.9676, 126.4514, 198.0641, 324.5155, 506.4506, 830.3209, 1185.804, 1840.641, 2322.576, 2799.994, 4116.121, 7129.018, 10322.56, 16774.16] }

    psvTable = pd.DataFrame(data, columns = ["Designation", "typicalFlanges", "areaIN2", "areaMM2"])


    df_mask = psvTable['areaMM2']>=myAreaMM2
    maskTable = psvTable[df_mask]  # now we can find the n smallest of these
    myPSV = maskTable.nsmallest(1,'areaMM2')
    return(myPSV.iloc[0]['Designation'])


def thermExpansionRate(heat, alpha, heatCap):
    # liquid thermal expansion sizing
    # heat kJ/s,
    # alpha thermal expansion coefficient, 1/K
    # heatCap kJ/kg.K
    # return mass flow(kg/hr) = 3600 * alpha (1/K) * heat (kJ/s) / (heatCap (kJ/kg.K))
    return (3600.0 * alpha * heat / heatCap)

def liquidVaporizeReliefRate(heat,latent):
    # heat kJ/s
    # latent is in kJ/kg
    # flow rate kg/h = (heat/1000) * (1/latent) * (1/3600)
    return (heat * (1.0/latent)* 3600.0)

def poolFireReliefRate(wettedAreaM2,latent,prompt):
    # wettedAreaM2 is in m2
    # latent is in kJ/kg
    # prompt fire response is either "prompt" or something else
    # heat Watts = C1 F A^0.82
    # flow rate kg/h = (heat/1000) * (1/latent) * (1/3600)

    C1 = 70900.0   # prompt fire fighting DOES NOT exist

    if (prompt == "prompt"):
        C1 = 43200.0   # prompt fire fighting exists
    
    F = 1.0 # no credits
    Q = C1*F*wettedAreaM2**.82
    flowRate = (Q/1000) * (1/latent)* 3600
    return (flowRate)

# steam flux needs a field for State: Saturated or a temperature in C
def PSVsteamFlux(Pkpa, State):
    # Napier equation for steam flow through PSV
    # mass flow = area*flux
    # flux will be in kg/(mm2.hr)
    # Good for choked flow of steam through orifice
    # Wlbhr = Ain2 Cnapier Kd Ppsi Ksh Kb Kn
    # Wkghr = (1/2.205)*(1/25.4^2)*Amm2 Cnapier Kd Ppsi Ksh Kb Kn
    # flux = (1/2.205)*(1/25.4^2) Cnapier Kd Ppsi Ksh Kb Kn
    # our function will return the flow in kg/h
    # is the napier constant 51.45 or something else?
    # must agree on the value for Kd, 0.975 or 0.938
    # discharge coefficient is usually 0.975 for steam and vapour but
    # could be different depending on manufacturer
    # flowRate = Area * Flux    
    Cnapier = 51.45
    Kd = 0.975
    Ksh = getKsh(Pkpa, State)
    Kb = 1.0
    Kn = 1.0     # low pressure and high pressure correction
    if (Pkpa > 10300.0):
        Kn = (2.7644*Pkpa/100.0 - 1000.0)/(3.3242*Pkpa/100.0 - 1061.0)
    
    Ppsi = Pkpa * (14.503773800721813/100)
    return (Cnapier*Kd*Ppsi*Ksh*Kb*Kn/(2.205*25.4**2))

def PSVsteamRate(areaMM2, Pkpa, State):
    # Napier equation for steam flow through PSV
    # Good for choked flow of steam through orifice
    # Pkpa is in kPaa
    # our function will return the flow rate in kg/hr
    # flowRate = Area * Flux
    psvFlux = PSVsteamFlux(Pkpa, State) # get the flux through the PSV, kg/hr.mm2
    Wkg = areaMM2*psvFlux
    return (Wkg)

def PSVsteamSize(Wkg, Pkpa, State):
    # Napier equation for steam flow through PSV
    # Good for choked flow of steam through orifice
    # our function will return the area in MM2
    # flowRate = Area * Flux
    psvFlux = PSVsteamFlux(Pkpa, State) # get the flux through the PSV, kg/hr.mm2
    Amm2 = Wkg/psvFlux
    return (Amm2)


def apiC(k):
    # SI form
    C = 0.03948*math.sqrt(k*(2.0/(k+1.0))**((k+1.0)/(k-1.0)))   # C coefficient API 520A fig 32
    return (C)


def PSVvaporFlux(P, Tcelcius, MW, k, Z):
    # vapour equation for PSV rating
    # P in kPaa
    # Tcelcius in celcius,
    # k is ratio of specific heats
    # Z compressibility factor
    # MW mole weight
    # Patm is atm pressure in kPaa
    # return value is flux in kg/mm2.hr
#    coeffSI = 13160.0; choop this out, not in the new version
    T = Tcelcius + 273.15;
    Kd = 0.975 # this could be reduced
    Kb = 1.0 # do not consider backpressure derating
    Kc = 1.0 # no derating for rupture disc
    C = apiC(k)
    rootTerm = math.sqrt(T*Z/MW);
    coeffTerm = 1.0 / (C * Kd * P * Kb * Kc)
    #    areaMM2 = W*coeffTerm*rootTerm;
    # W = areaMM2 * flux
    # flux = 1.0/(CoeffTerm * rootTerm)

    return (1.0/(coeffTerm*rootTerm))


def PSVvaporRate(areaMM2, P, Tcelcius, MW, k, Z):
    # vapour equation for PSV rating
    # areaMM2 orifice area mm2
    # W flow rate in kg/h
    # P in kPaa
    # Tcelcius in celcius,
    # k is ratio of specific heats
    # Z compressibility factor
    # MW mole weight
    # Patm is atm pressure in kPaa
    # return value is rated flow rate kg/h
    # W = areaMM2 * flux
    # flux = 1.0/(CoeffTerm * rootTerm)
    psvFlux = PSVvaporFlux(P, Tcelcius, MW, k, Z) # get the flux through the PSV, kg/hr.mm2

    return (areaMM2 * psvFlux)

def PSVvaporSize(W, P, Tcelcius, MW, k, Z):
    # vapour equation for PSV sizing
    # W flow rate in kg/h
    # P in kPaa
    # Tcelcius in celcius,
    # k is ratio of specific heats
    # Z compressibility factor
    # MW mole weight
    # Patm is atm pressure in kPaa
    # return value is required area mm2
    # W = areaMM2 * flux
    # flux = 1.0/(CoeffTerm * rootTerm)
    psvFlux = PSVvaporFlux(P, Tcelcius, MW, k, Z) # get the flux through the PSV, kg/hr.mm2

    return (W / psvFlux)




def PSVliquidSize(W, P, Pback, d, mu):
    # liquid equation for PSV sizing, PSV requiring capacity certification
    # W flow rate in kg/h
    # P in kPag
    # Pback in kPag
    # d density kg/m3
    # mu viscosity cP or mPa.s
    # Patm is atm pressure in kPaa
    # return value is required area mm2
    convergeEPS = 1.0e-4 # converge the viscous correction to this value
    coeffSI = 11.78;
    Q = 1000.0*W/(d*60.0);  # litres/minute
    G = d/1000.0; # specific gravity
    P1 = P;
    P2 = Pback;
    
    Kd = 0.65 # This is for a PSV. If we look for a rupture disk, Kd = 0.62
    Kw = 1.0 # do not consider backpressure derating
    Kc = 1.0 # no derating for rupture disc
    
    Kv = 1.0 # first pass at viscous correction, this will be iterated
    i = 1
    areaMM2 = 0.0
    maxIter = 10
    errConverge = 1.0 
    maxIter = 10
    while (errConverge > convergeEPS) and ( i < maxIter):   
        rootTerm = math.sqrt(G/(P1-P2));
        coeffTerm = coeffSI * Q / (Kd * Kw * Kc * Kv)
        areaMM2 = coeffTerm*rootTerm;
#        diam = sqrt(4.0*areaMM2/pi)/1000;
        R = Q * 18800 * G / (mu*math.sqrt(areaMM2)); # Reynolds number
        oldKv = Kv
        Kv = 1.0 / (0.9935 + 2.878/math.sqrt(R) + 342.75/(R**1.5))
        errConverge = abs(oldKv - Kv)

        i += 1
#    return (i)  # mass flow rate kg/h
    return (areaMM2)

def PSVliquidRate(areaMM2, P, Pback, d, mu):
    # liquid equation for PSV flow rate
    # W flow rate in kg/h
    # areaMM2 PSV area
    # P in kPag
    # Pback in kPag
    # d density kg/m3
    # mu viscosity cP or mPa.s
    # Patm is atm pressure in kPaa
    # return value is flow rate kg/h
    # I have not edited this function yet
    # from API 520
    # A = (11.78 Q) / (Kd Kw Kc Kv) * sqrt(G/DP)
    # Q = (A Kd Kw Kc Kv/11.78) sqqrt(DP/G)
    coeffSI = 11.78;
#    Q = 1000*W/(d*60.0);  # litres/minute
#    A = areaMM2 / (25.4^2)
    convergeEPS = 1.0e-4 # converge the viscous correction to this value
    G = d/1000; # specific gravity
    P1 = P;
    P2 = Pback;

    Kd = 0.65 # typical value
    Kw = 1.0 # do not consider backpressure derating
    Kc = 1.0 # no derating for rupture disc

    Kv = 1.0 # first pass at viscous correction, this will be iterated
    i = 1
    Q = 0.0 # initialize the variable
    W = 0.0;  # initialize the variable
    maxIter = 10
    errConverge = 1.0 
    while (errConverge > convergeEPS) and ( i < maxIter):   
        # areamm2 = coeff * root
        rootTerm = np.sqrt((P1-P2)/G);
        coeffTerm = (Kd * Kw * Kc * Kv) / coeffSI
        Q = areaMM2*coeffTerm*rootTerm;   # litres per minute
#        diam = sqrt(4.0*areaMM2/pi)/1000.0;
        R = Q * 18800.0 * G / (mu*math.sqrt(areaMM2));  # Reynolds number
        oldKv = Kv
        Kv = 1.0 / (0.9935 + 2.878/math.sqrt(R) + 342.75/(R**1.5))
        errConverge = abs(oldKv - Kv)
        i += 1
        
    return (Q*60*d/1000.0)  # mass flow rate kg/h
#    return (i)  # mass flow rate kg/h




