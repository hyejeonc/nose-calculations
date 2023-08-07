% Author: Oivind Wilhelmsen, Date: 2012-11-21
% This function calculates the ideal gas chemical potential of the 
% components with the derivatives, NB, does not take arrays. 

function mu = IG_mu_a(T,P,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units      Dim.     Type
%
% INPUT:   T        Temperature                K        1 x 1    double
%          P        Pressure                   Pa       1 x 1    double
%          Z        Mole fractions             (-)      1 x nc   double
%
% OUTPUT:  mu       Chemical potential         J/kmol   n x nc   struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

if length(T)>1
    disp('Please send in just on state')
    stop
end

% The pure component ideal gas quantities [n x nc]
S0 = IG_S0_a(T);         % Entropies
H0 = IG_H0_a(T);         % Enthalpies

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

% Chemical potential of the mixture [J/ kmole]
mu.mu = H0-T.*S0 + ktp1.R.*T.*(log(P./ktp1.Pref)+term);

% Derivatives
mu.dT = -1*S0+ktp1.R.*(log(P./ktp1.Pref)+term);  % Temperature [J/kmole K]
mu.dP = T*ktp1.R./P;                            % Pressure [J/kmole P]

% Scaled chemical potential derivative n_tot*(dmu_ig/dni) [J/kmole]
for i=1:ktp1.nc
    for j=1:ktp1.nc
        
        if i==j
        mu.dZ(i,j) = ktp1.R.*T.*(1./Z(i)).*(1-Z(i));
        else
        mu.dZ(i,j) = ktp1.R.*T.*(1./Z(i)).*(-Z(i));
        end
        
    end
end