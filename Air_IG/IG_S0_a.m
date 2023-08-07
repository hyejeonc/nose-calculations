% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas entropy of each 
% components using a fourth degree polynomial in temperature

function S0 = IG_S0_a(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units      Dim.     Type
%
% INPUT:   T        Temperature                K        n x 1    double
%
% OUTPUT:  S0       Ideal gas entropy          J/kmol K nc x 1   double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

if min(T)<max(ktp1.Cp_Tmin) || max(T)>min(ktp1.Cp_Tmax)
    if ktp1.Flag_print
    disp(' NB! The ideal gas CP-polynomial is not valid')
    end
end

% The difference between the temperature and the reference [K]

for i=1:ktp1.nc
    
    % Compute the integrated Cp/T;
    S0(i,:) = ktp1.Cp(i,1)*(log(T)-log(ktp1.Tref)) + ... 
              ktp1.Cp(i,2)*(T-ktp1.Tref)           + ...
             (ktp1.Cp(i,3)/2)*(T.^2 - ktp1.Tref^2) + ...
             (ktp1.Cp(i,4)/3)*(T.^3 - ktp1.Tref^3) + ...
             (ktp1.Cp(i,5)/4)*(T.^4 - ktp1.Tref^4);
    
    % Include the enthalpy standard state     
    S0(i,:)=S0(i,:)*ktp1.R+ktp1.S_std(i);     
end