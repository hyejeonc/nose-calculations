% Author: Oivind Wilhelmsen, Date: 2012-12-19
% Calculates the chemical potential from cubic EOS given temperature,
% pressure and composition and derivatives 
% with respect to temperature, concentration and composition:

function [mu_struct] = TP_ChemicalPotential_a(T,P,Z,flag_deriv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol     Explanation              Units    Dim.     Type
%
% INPUT:   T          Temperature                K      1 x 1    double
%          P          Pressure                   Pa     1 x 1    double
%          Z            Mole fractions             -        -      struct
%          Phase        Phase flag, 1=liq, 2=gas   -      1 x 1      -%          
%          flag_deriv   Derivative flag:           -        -      integer 
%                       0=no derivatives
%                       1=all derivatives
%
% OUTPUT:  mu_struct  Chemical potential         -        -      struct
%                     struct     
%          Temp       Struct with temporary      -        -      struct
%                     variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

phase   = 2;

% Product [J/kmole]
Prod=ktp1.R*T;

% Calculate the ideal gas contribution
mu_struct_ig = IG_mu_a(T,P,Z);

% Calculate the residual contribution ln(Phi(i))=mu_res/RT
Fug_struct =TP_Fug_a(T,P,Z,phase,flag_deriv);

% The chemical potential [J/kmol]
mu_struct.mu=(mu_struct_ig.mu+Fug_struct.Ln_Fug*Prod)./ktp1.Mw;

if flag_deriv>0
    % The temperature derivative [J/kmol K]
    mu_struct.dT=(mu_struct_ig.dT+Fug_struct.Ln_Fug*ktp1.R+...
             Fug_struct.dT*Prod)./ktp1.Mw;
            
    % The concentration derivative [J/m^3]
    mu_struct.dP=mu_struct_ig.dP+Fug_struct.dP*Prod;

    % The scaled composition derivative [J/kmol]
    mu_struct.dN=mu_struct_ig.dZ+Fug_struct.dN*Prod;
end
