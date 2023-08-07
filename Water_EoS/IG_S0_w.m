% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas entropy of each 
% components using a fourth degree polynomial in temperature

function S0 = IG_S0_w(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units      Dim.     Type
%
% INPUT:   T        Temperature                K        n x 1    double
%
% OUTPUT:  S0       Ideal gas entropy          J/kmol K nc x 1   double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp2

if min(T)<max(ktp2.Cp_Tmin) || max(T)>min(ktp2.Cp_Tmax)
    if ktp2.Flag_print
    disp(' NB! The ideal gas CP-polynomial is not valid')
    end
end

% The difference between the temperature and the reference [K]

for i=1:ktp2.nc
    
    % Compute the integrated Cp/T;
    S0(i,:) = ktp2.Cp(i,1)*(log(T)-log(ktp2.Tref)) + ... 
              ktp2.Cp(i,2)*(T-ktp2.Tref)           + ...
             (ktp2.Cp(i,3)/2)*(T.^2 - ktp2.Tref^2) + ...
             (ktp2.Cp(i,4)/3)*(T.^3 - ktp2.Tref^3) + ...
             (ktp2.Cp(i,5)/4)*(T.^4 - ktp2.Tref^4);
    
    % Include the enthalpy standard state     
    S0(i,:)=S0(i,:)*ktp2.R+ktp2.S_std(i);     
end