% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas heat capacity of 
% the components using a fourth degree polynomial in temperature

function Cp0 = IG_Cp0_a(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      n x 1    double
%
% OUTPUT:  Cp0      Heat capacity           J/kmol K  nc x 1   double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

if min(T)<max(ktp1.Cp_Tmin) || max(T)>min(ktp1.Cp_Tmax) 
    if ktp1.Flag_print
    disp(' NB! The ideal gas CP-polynomial is not valid')
    end
end

for i=1:ktp1.nc
    Cp0(i,1) = ktp1.Cp(i,1) + ktp1.Cp(i,2)*T + ktp1.Cp(i,3)*T.^2 + ...
              ktp1.Cp(i,4)*T.^3 + ktp1.Cp(i,5)*T.^4; 
end

% Use the gas constant to obtain correct units [J/kmol K]
Cp0=Cp0*ktp1.R;