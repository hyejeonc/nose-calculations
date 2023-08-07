% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas entropy of the mixture and
% the derivatives using a fourth degree polynomial in temperature

function S = IG_S(T,P,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units      Dim.     Type
%
% INPUT:   T        Temperature                K        n x 1    double
%          P        Pressure                   Pa       n x 1    double
%          Z        Composition                (-)      n x nc   double
%
% OUTPUT:  S        Ideal gas entropy          J/kmol K ncx1     struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

% The entropies of the components [n x nc]
S0 = IG_S0_a(T);

% Calculate the mixing term
for i = 1:length(Z(:,1))
  for j = 1:length(Z(1,:))
    if Z(i,j)>10^(-15)
      term(i,j) = log(Z(i,j));
    else
      term(i,j) = 0;
    end
  end
end

% Entropy of the mixture [J/K  mole]
S.S = sum(Z .* S0) -ktp1.R * sum(Z.*term) -ktp1.R.*log(P./ktp.Pref);

% Derivatives
S.dT=IG_Cp_a(T,Z)./T;                           % Temperature derivatives [J/kmol K^2]
S.dP=-ktp1.R./P;                               % Pressure derivative [J/kmol K Pa]

for i=1:ktp.nc
    S.dZ(i,:)=S0(i,1)-ktp.R.*term(i,1)-ktp1.R.*log(P./ktp1.Pref);  % Composition derivative [J/kmol K]
end