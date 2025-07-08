function [DIC,HCO3,CO3,pH,pCO2] = f_csys_alk_CO2(TC,S,ALK,CO2)
%
% 	Function based on files: csys.m and equic.m   
%   by Richard E. Zeebe & Dieter A. Wolf-Gladrow
%
%   Modified by Katja Fennel to work as a function.
%
%   Using their flag 4 (alk and CO2 given)
%
%   Purpose: Given two components of the CO2 system, calculate all other.
%   In this case, alk and CO2 are given at pressure zero (surface ocean).
%
%   Input parameters:
%         TC :  temperature (degrees C); valid range: ~0-35 deg C
%         S :   salinity; valid range: ~20-40
%         ALK : alkalinity (mumol equivalents/kg)
%         CO2 : aqueous CO2 (mumol/kg)
%
%   Output parameters:
%         DIC : dissolved inorganic carbon (mumol/kg)
%         HCO3 : bicarbonate ion concentration (mumol/kg)
%         CO3 : carbonate ion concentration (mumol/kg)
%         pH :  pH (total scale)
%         pco2 : partical pressure of CO2 (ppm)
%
% ---------------------------------------------------------------
%----  define constants
T = TC + 273.15;        % TC [C]; T[K]; conversion [deg C] <-> [K])
R = 83.14510;           % mol bar deg-1; conversion fco2 - pco2

% --------------------- Kwater -----------------------------------
%       Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
%       K_w in mol/kg-soln.
%       pH-scale: pH_total ('total` scale).
%                                                       
lnKw = -13847.26/T + 148.9652 - 23.6521*log(T) ...
      + (118.67/T - 5.977 + 1.0495*log(T))*sqrt(S) - 0.01615*S;
Kw  = exp(lnKw);

%---------------------- Kh (K Henry) ----------------------------
%               CO2(g) <-> CO2(aq.)
%               Kh      = [CO2]/ p CO2
%
%   Weiss (1974)   [mol/kg/atm]
%
%   Weiss + Murray & Riley (1971) data: 
%   T = 1-35 deg C, S = 0-38.
%                             
%
nKhwe74 = 9345.17 / T - 60.2409 + 23.3585 * log(T/100.) ...
          + S*(0.023517-0.00023656*T+0.0047036e-4*T*T);
Kh= exp(nKhwe74);

%---------------------- p2f (pCO2 -> fCO2) -----------------------               
%   convert from pCO2 to fCO2
%   Weiss (1974)
%
Pstd = 1.01325;
delC = (57.7 - 0.118*T);
B = -1636.75 + 12.0408*T - 0.0327957*T^2 + 3.16528*0.00001*T^3;
p2f = exp((B + 2*delC)*Pstd/(R*T));

% --------------------- K1 ---------------------------------------
%   first acidity constant:
%   [H^+] [HCO_3^-] / [H_2CO_3] = K_1
%
%   Mehrbach et al (1973) refit by Lueker et al. (2000).
%
%   Mehrbach data: T = 2-35 deg C, S = 19-43.
%
%   pH-scale: 'total'. mol/kg-soln
pK1 = 3633.86/T - 61.2172 + 9.6777*log(T) - 0.011555*S + 0.0001152*S*S;
K1  = 10^(-pK1);

% --------------------- K2 ----------------------------------------
%   second acidity constant:
%   [H^+] [CO_3^--] / [HCO_3^-] = K_2
%
%   Mehrbach et al. (1973) refit by Lueker et al. (2000).
%
%   Mehrbach data: T = 2-35 deg C, S = 19-43.
%
%   pH-scale: 'total'. mol/kg-soln
pK2 = 471.78/T + 25.9290 - 3.16967*log(T) - 0.01781*S + 0.0001122*S*S;
K2  = 10^(-pK2);

% --------------------- Kb  --------------------------------------------
%
%  Kbor = [H+][B(OH)4-]/[B(OH)3]
%
%   (Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
%   pH-scale: 'total'. mol/kg-soln
tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S^(3/2.)-0.0996*S*S);
tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S)*log(T);
Kb = exp(tmp1 / T + tmp2 + tmp3 + 0.053105*sqrt(S)*T);

bor = 1*(416*(S/35.))*1.e-6;   % (mol/kg), DOE94

%
% ------  actual calculations
s = CO2*1.e-6;
alk = ALK*1.e-6;

p4 = 1.;              
p3 = Kb+alk;
p2 = alk*Kb-s*K1-Kb*bor-Kw;
p1 = -s*Kb*K1-s*2.*K1*K2-Kw*Kb;
p0 = -2.*s*Kb*K1*K2;
p = [p4 p3 p2 p1 p0];
r = roots(p);
h = max(real(r));

dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
fco2 = s/Kh;
pco2 = fco2/p2f;

% ----------- change units: mumol/kg
pH      = -log10(h);
pCO2    = pco2*1.e6;
HCO3    = hco3*1.e6;
CO3     = co3*1.e6;
DIC     = dic*1.e6;

return;

