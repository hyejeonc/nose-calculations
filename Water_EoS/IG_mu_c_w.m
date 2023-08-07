% Author: Oivind Wilhelmsen, Date: 2012-11-27
% This function calculates the ideal gas chemical potential of the 
% components with the derivatives, NB, does not take arrays. 
% This version takes the molar concentration as input variable

function mu = IG_mu_c(T,c,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units      Dim.     Type
%
% INPUT:   T        Temperature                K        1 x 1    double
%          c        Concentration              mol/m^3  1 x 1    double
%          Z        Mole fractions             (-)      nc x 1   double
%
% OUTPUT:  mu       Chemical potential         J/kmol   n x nc   struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp2

% vector to create the right dimension:
unit_vec=ones(size(Z));

if length(T)>1
    disp('Please send in just on temperature and pressure')
    stop
end

% The pure component ideal gas quantities [n x nc]
S0 = IG_S0_w(T);         % Entropies
H0 = IG_H0_w(T);         % Enthalpies

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

% Calculate the pressure [Pa]
P=c*ktp2.R*T;

% Chemical potential of the mixture [J/ kmole]
mu.mu = H0-T.*S0 + ktp2.R.*T.*(log(P./ktp2.Pref)+term);

% Derivatives (constant c)
mu.dT = -1*S0+ktp2.R.*...
        (log(P./ktp2.Pref)+term)+ktp2.R;      % Temperature [J/kmole K]
mu.dc = unit_vec*T*ktp2.R./c;                % Composition [J/kmole P]
mu.dV = -1*unit_vec*T*ktp2.R.*c;             % Volume [J/m^3]

% Scaled chemical potential derivative n_tot*(dmu_ig/dni) (given c) [J/kmole]
for i=1:ktp2.nc
    for j=1:ktp2.nc
        
        if i==j
        mu.dN(i,j) = ktp2.R.*T.*(1./Z(i)).*(1-Z(i));
        else
        mu.dN(i,j) = ktp2.R.*T.*(1./Z(i)).*(-Z(i));
        end
        
    end
end

% Scaled chemical potential derivative n_tot*(dmu_ig/dni) (given V) [J/kmole]
mu.dN_v=mu.dN+ktp2.R*T;