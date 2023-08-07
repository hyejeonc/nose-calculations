% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas enthalpy of each
% components using a fourth degree polynomial in temperature

function H0 = IG_H0_a(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      n x 1    double
%
% OUTPUT:  H0       Ideal gas enthalpy        J/kmol  nc x 1   double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

if min(T)<max(ktp1.Cp_Tmin) || max(T)>min(ktp1.Cp_Tmax)
    if ktp1.Flag_print
    disp(' NB! The ideal gas CP-polynomial is not valid')
    end
end

for i=1:ktp1.nc
    
    % Compute the integrated polynomial
    H0(i,:) = ktp1.Cp(i,1)*(T-ktp1.Tref) + ...
              (ktp1.Cp(i,2)/2)*(T.^2-ktp1.Tref.^2) +...
              (ktp1.Cp(i,3)/3)*(T.^3-ktp1.Tref.^3) + ...
              (ktp1.Cp(i,4)/4)*(T.^4-ktp1.Tref.^4) + ...
              (ktp1.Cp(i,5)/5)*(T.^5-ktp1.Tref.^5);
    
    % Include the enthalpy standard state     
    H0(i,:)=H0(i,:)*ktp1.R+ktp1.H_std(i);     
end

