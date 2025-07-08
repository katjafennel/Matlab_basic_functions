# README

This folder contains functions for calculating seawater carbonate system parameters. I recommend adding these functions to a folder in your Matlab search path or, even better, automatically adding this folder to your search path whenever you start up Matlab by including this line
```
addpath /insert-path-to-this-folder-on-your_system
``` 
in your `startup.m` file.

## Modified csys functions (Zeebe and Wolf-Gladrow)

The following functions are based on scripts `csys.m` and `equic.m` provided by Richard E. Zeebe and Dieter A. Wolf-Gladrow as code supplement for the textbook CO2 in seawater: equilibrium, kinetics, isotopes (Elsevier Oceanographic Series, Vol. 65. 346 pp, 2001).
  
```
f_csys_alk_CO2.m
f_csys_alk_DIC.m
f_csys_alk_pCO2.m
f_csys_h_CO2.m
f_csys_h_pCO2.m
```

The scripts by Zeebe and Wolf-Gladrow were modified by Katja Fennel to work as a functions. Given two components of the CO2 system, they calculate all other. For details on input and output parameters, use the built-in help function, i.e. call `>> help f_csys_alk_CO2.m`.

## CO2 solubility and Schmidt number

The following functions calculate the solubility of CO2 and the Schmidt number using coefficients from R.F. Weiss (Marine Chemistry 2(3), 203-215, 1974). These are needed for calculating air-sea flux of CO2.
  
```
Ko_Weiss.m
SCO2_Weiss.m
```

For details on input and output parameters, use the built-in help function.


## History
The content of the README file is up to date as of July 7, 2025. -KF





