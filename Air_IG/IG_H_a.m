% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas enthalpy of the
% mixture using a fourth degree polynomial in temperature

function H = IG_H_a(T,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      n x 1    double
%          Z        Compositions               -      n x nc   double
%
% OUTPUT:  H        Ideal gas enthalpy        J/kmol  nc x n   struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The enthalpies of the components [nc x n]
H0 = IG_H0_a(T);

% Enthalpy of the mixture 
H.H=sum(Z.*H0);

% Derivatives
H.dT=IG_Cp_a(T,Z); % Temperature derivatives [J/kmol K]
H.dP=0;          % Pressure derivative [J/kmol Pa]
H.dZ=H0;         % Composition derivative [J/kmol]