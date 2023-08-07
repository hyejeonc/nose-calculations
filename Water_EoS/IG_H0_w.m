% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas enthalpy of each
% components using a fourth degree polynomial in temperature

function H0 = IG_H0_w(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      n x 1    double
%
% OUTPUT:  H0       Ideal gas enthalpy        J/kmol  nc x 1   double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp2

if min(T)<max(ktp2.Cp_Tmin) || max(T)>min(ktp2.Cp_Tmax)
    if ktp2.Flag_print
    disp(' NB! The ideal gas CP-polynomial is not valid')
    end
end

for i=1:ktp2.nc
    
    % Compute the integrated polynomial
    H0(i,:) = ktp2.Cp(i,1)*(T-ktp2.Tref) + ...
              (ktp2.Cp(i,2)/2)*(T.^2-ktp2.Tref.^2) +...
              (ktp2.Cp(i,3)/3)*(T.^3-ktp2.Tref.^3) + ...
              (ktp2.Cp(i,4)/4)*(T.^4-ktp2.Tref.^4) + ...
              (ktp2.Cp(i,5)/5)*(T.^5-ktp2.Tref.^5);
    
    % Include the enthalpy standard state     
    H0(i,:)=H0(i,:)*ktp2.R+ktp2.H_std(i);     
end

