% Author: Oivind Wilhelmsen, Date: 2012-11-28
% Calculates the chemical potential from cubic EOS given temperature,
% total concentration and composition and derivatives 
% with respect to temperature, concentration and composition:

function mu_struct  = TC_ChemicalPotential_a(T,c,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol     Explanation              Units    Dim.     Type
%
% INPUT:   T          Temperature                K      1 x 1    double
%          c          Molar concentration        kmol/m3 1 x 1    double
%          Z          Mole fractions             -    noc x 1    double
%
% OUTPUT:  mu_struct  Chemical potential         -        -      struct
%                     struct     
%          Temp       Struct with temporary      -        -      struct
%                     variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

% Product [J/kmole]
Prod=ktp1.R*T;

% Calculate the ideal gas contribution
mu_struct_ig = IG_mu_c_a(T,c,Z);

% Calculate the residual contribution to the Helmholtz energy
[~, F_struct] =TC_Ares_a(T,c,Z);

% The chemical potential [J/kmol]
mu_struct.mu=mu_struct_ig.mu+F_struct.mures*Prod;

% The temperature derivative [J/kmol K]
mu_struct.dT=mu_struct_ig.dT+F_struct.mures*ktp1.R+...
             F_struct.dT*Prod;
         
% The concentration derivative [J m^3/kmol^2]
mu_struct.dc=mu_struct_ig.dc+F_struct.dc*Prod;

% The volume derivative [J/kmol m^3]
mu_struct.dV=mu_struct_ig.dV+F_struct.dV*Prod;

% The scaled composition derivative [J/kmol]
mu_struct.dN=mu_struct_ig.dN+F_struct.dN*Prod;

% The scaled composition derivative constant volume [J/kmol]
mu_struct.dN_v=mu_struct_ig.dN_v+F_struct.dN_v*Prod;