function Cp0 = IG_Cp_w(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      n x 1    double
%
% OUTPUT:  Cp0      Heat capacity           J/kmol K  nc x 1   double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp2

if min(T)<max(ktp2.Cp_Tmin) || max(T)>min(ktp2.Cp_Tmax) 
    if ktp2.Flag_print
    disp(' NB! The ideal gas CP-polynomial is not valid')
    end
end

for i=1:ktp2.nc
    Cp0(i,1) = ktp2.Cp(i,1) + ktp2.Cp(i,2)*T + ktp2.Cp(i,3)*T.^2 + ...
              ktp2.Cp(i,4)*T.^3 + ktp2.Cp(i,5)*T.^4; 
end

% Use the gas constant to obtain correct units [J/kmol K]
Cp0=Cp0*ktp2.R;