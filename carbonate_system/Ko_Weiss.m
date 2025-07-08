function [Ko]=Ko_Weiss(T,S)
%
% Solubility of CO2 in Seawater
% from Weiss, R.F., Marine Chemistry 2(3), 203-215 (1974)
%
% Inputs:
%         T [degC] temperature in Celsius
%         S [] salinity
%
% Output:
%         Ko [mol /kg /atm] solubility 
%
A = [-60.2409, 9345.17, 23.3585];       % mol /kg /atm
B = [0.023517, -0.00023656, 0.0047036]; % mol /kg /atm
T = T+273.15;                           % convert temperature to Kelvin
ln_Ko = A(1)+A(2)/T+(A(3)*log(T/100))+S*(B(1)+B(2)*T+B(3)*(T/100).^2);
Ko = exp(ln_Ko);
end