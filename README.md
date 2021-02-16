<!-- ABOUT THE PROJECT -->

## About The Project

This project provides functions in Python for standard calculations for Pressure Safety Valve (PSV) flow rate and sizing. The intention is to use the routines in a Jupyter Notebook file for documenting engineering work.  

The calculations should be adequate for engineering consulting work and preliminary sizing or rating. Definative sizing or rating calculations should be performed with methodologies or rating factors provided by the PSV manufacturer.


### Built With

The code is written in Python (and is the authors first Python package uploaded to GitHub). The code is intended to be used in a Jupyter Notebook. I have not used the routines in a stand-alone Python environment.


<!-- GETTING STARTED -->
## Getting Started

The following lines of code are needed in a Jupyter Notebook (Python shell) to pull the _unregistered_ package from GitHub and use the package.

~~~~
TBD
~~~~

### Prerequisites

The package requires the following packages

* numpy
* math
* scipy

<!-- TESTING -->
### Testing

The following code tests are available

* steam: size a steam PSV
* vapour: size a vapour PSV
* liquid1: size a liquid PSV
* liquid2: rate the capacity of a liquid PSV


<!-- USAGE EXAMPLES -->
## Usage

Refer to the Jupyter Notebook file for an example of how the code is used.

_For more examples, please refer to the [Documentation](https://example.com)_

Functions are provided for rating and sizing of Steam, Vapour and Liquid PSVs. Units of measure used in the package are:

* Flow rate, kg/h
* PSV flow area, mm2
* Pressure, kPaa
* Temperature, deg C

Functions

* PSVsteamRate(areaMM2, Pkpa, State)
    * areaMM2 is the API flow area in mm2
    * State is either a temperature in deg C (superheated steam) or a string ("Sat", saturated steam).
    * return value is kg/h
* PSVsteamSize(Wkg, Pkpa, State)
    * Wkg is mass flow rate in kg/h
    * State is either a temperature in deg C (superheated steam) or a string ("Sat", saturated steam).
    * return value is orifice area mm2
* PSVsteamFlux(Pkpa, State)
    * this is the main function for steam PSV calculations
    * the return value is mass flux in kg/hr.mm2
    * this is used to calculate either the area (given the flow rate) or the flow rate (given the area)

    $$
    K_d = 0.975 \\
	K_{sh} = \mbox{Superheat derating (lookup table)} \\
    K_b = 1.0 \, \mbox{(no backpressure derating)} \\
    K_n = \frac{2.7644 \times Pkpa/100.0 - 1000.0}{3.3242 \times Pkpa/100.0 - 1061.0}, P > 10300 \mbox{kPa} \\
    Ppsi = Pkpa \times (14.503773800721813/100) \\
    flux_{kg/hr.mm^2} = 51.45 \times K_d \times Ppsi \times K_{sh} \times K_b \times K_n / (2.205 \times 25.4^2)
    $$ 
    
* PSVvaporRate(areaMM2, P, Tcelcius, MW, k, Z)
    * given the API orifice area, pressure (kPaa), temperature (deg C), mole weight, ratio of specific heats or isentropic coefficient, and compressibility factor
    * return value is the flow rate in kg/h
* PSVvaporSize(W, P, Tcelcius, MW, k, Z)
    * given the flow rate in kg/h, and the other standard inputs
    * return value is the PSV flow area in mm2
* PSVvaporFlux(P, Tcelcius, MW, k, Z)
    * this is the main function for vapour PSV calculations
    * the return value is mass flux in kg/hr.mm2
    * this is used to calculate either the area (given the flow rate) or the flow rate (given the area)

$$
    Kd = 0.975 \, \mbox{discharge coefficient, can vary with mfg} \\
    Kb = 1.0 \, \mbox{do not consider backpressure derating} \\
    Kc = 1.0 \, \mbox{no derating for rupture disc} \\
    C = 0.03948 \sqrt{  k \left(\frac{2.0}{k+1}\right)^{(k+1)/(k-1)}    }   \, \mbox{API 520A fig 32} \\
    flux_{kg/hr.mm2} = \frac{C * Kd * P * Kb * Kc}{\sqrt{T_{kelvin} \times Z/MW}}
$$

* PSVliquidRate(areaMM2, P, Pback, d, mu)
    * given the API orifice area in mm2, inlet pressure kPag, backpressure kPag, density in kg/m3 and viscosity in mPa.s (cP)
    * the return value is the liquid flow rate in kg/h

    $$
    Kd = 0.65 \\
    Kw = 1.0 \\
    Kc = 1.0 \\
        Q = A_{mm2} \frac{Kd Kw Kc Kv}{11.78} \sqrt{\Delta P/(\rho/1000)} \, \mbox{litres per minute} \\
        \mbox{where} \\
        R = \frac{Q * 18800 * \rho/1000}{\mu*\sqrt{A_{mm2}}} \, \\
        Kv = 1.0 / (0.9935 + 2.878/\sqrt{R} + 342.75/(R^{1.5})) \\
    flowrate = Q*60*\rho/1000.0 \,  \mbox{mass flow rate kg/h}
    $$
    
* PSVliquidSize(W, P, Pback, d, mu)
    * given the liquid relief rate in kg/h
    * return value is the PSV flow area in mm2

    $$
    Q = 1000.0*W/(\rho*60.0) \, \mbox{l/min} \\
    Kd = 0.65  \\
    Kw = 1.0 \\
    Kc = 1.0 \\
        A_{mm2} = 11.78 \frac{Q}{Kd * Kw * Kc * Kv} * \sqrt{(\rho/1000)/\Delta P} \\
        \mbox{where} \\
        R = \frac{Q * 18800 * \rho/1000}{\mu*\sqrt{A_{mm2}}} \, \\
        Kv = 1.0 / (0.9935 + 2.878/\sqrt{R} + 342.75/(R^{1.5})) 
    $$

Utility functions

* getKsh(PkPa, State)
    * this provides the superheat correction factor for the Napier steam formula
    * P: 140 - 20,600 kPa; T: saturated to 565 C
    * the pressure is in kPaa, and State is either a number (temperature in deg C) or a string (ie "Sat")
    * a simple lookup table is used to find the superheat value
* PSVareaOrifice(letter)
    * given the API orifice letter designation (D, E, F...)
    * return the API flow area in mm2
* PSVfindOrifice(area)
    * given the required API flow rate (from a sizing calculation)
    * return the letter designation for the next larger orifice size
* waterPsat(T_c)
    * given the temperature in deg C
    * return the saturation pressure of water in kPaa
    * this uses a hand correlated function of the form ln P_Pa = A + B/Tk + C*Tk + D*ln(Tk)
    * The data was fit in the pressure range 100 - 20,000 kPaa
    * a better approach would be to use IFC97 steam tables (future work)
* waterTsat(P_kPa)
    * given the pressure in kPaa
    * return the saturation temperature of water in deg C
    * this uses a hand correlated function of the form 1/T = A + B*ln P + C/ln P
    * The data was fit in the pressure range 100 - 20,000 kPaa
    * a better approach would be to use IFC97 steam tables (future work)
* thermExpansionRate(heat, alpha, heatCap)
    * calculate the relief flow rate for thermal expansion of a liquid
    * heat is the applied heat in kJ/s
    * alpha is the cubic thermal expansion coefficient, 1/K
    * heatCap is the heat capacity at constant pressure, kJ/kg.K
    * the return value is the mass flow rate in kg/h

    $$
    m = \frac{3600 * \alpha * q}{C_p} \, \mbox{kg/h}
    $$

* poolFireReliefRate(wettedAreaM2,latent,prompt)
    * this calculates the relief flow rate for a liquid vapourized in a fire scenario
    * wetted area refers to the wet surface area in the vessel, in m2
    * latent is the latent heat in kJ/kg
    * prompt is a string to denote fire response time: "prompt" or anything else

$$    
F = 1.0 \, \mbox{no environmental credits} \\
    flowRate = \frac{C1 * F * A_{w,m2}^{0.82}}{1000 \Delta H} * 3600 \, \mbox{kg/h} \\
    \mbox{where} \\
        C1 = 70900 \, \mbox{if prompt fire fighting DOES NOT exist} \\
        C1 = 43200 \, \mbox{if prompt fire fighting exists} 
    $$
    
* liquidVaporizeReliefRate(heat,latent)
    * this calculates the relief rate for vapourizing a liquid with a specified heat source
    * heat supply is in kJ/s

    $$
    m = \frac{q}{\Delta H} * 3600 \, \mbox{kg/h}
    $$

Refer to the Jupyter notebook file PSVreliefExample.ipynb for working examples.

<!-- ROADMAP -->
## Roadmap

* implement IAPWS IFC97 steam tables (ie Xsteam)
* implement Homogeneous Equilibrium Method (HEM) for steam. This permits relief calculations for two phase water-steam mixtures.
* provide ability to over-ride the standard PSV constants, such as discharge coefficient.
* implement backpressure correction for standard and balanced bellows PSV.



<!-- CONTRIBUTING -->
## Contributing

Send me a note.



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Kevin Dorma - [@kevindorma](https://twitter.com/KevinDorma) - kevin@kevindorma.ca

Project Link: [https://github.com/kevindorma/psvpy](https://github.com/kevindorma/psvpy)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

Not sure who to acknowledge.