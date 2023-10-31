% Author: Oivind Wilhelmsen, Date: 2012-11-22
% This function calculates the omega_A and omega_B parameters for a cubic 
% EoS using the requirement that (dP/dV)=0 and (d2P/dV^2)=0 at the
% critical point

function k = CB_constants(k, EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation                   Dimensions    Type
%
% INPUT:   k        Constants for the calculations     -        struct
%
% OUTPUT:  k        Constants for the calculations     -        struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum number of iterations
MaxIter=10;

% Convergence criteria
Eps=1.0E-10;

% Initiate and start third order method for Zc/omega_B
Zcomb=5.0;

    
% Start Newtons method for determination of consistent omega_A and omega_B

for i=1:MaxIter

    N  = Zcomb - 1.0;
    N1 = Zcomb - k.EoS.m1;
    N2 = Zcomb - k.EoS.m2;
    F  = N*N1^2+N*N1*N2+N*N2^2-N1^2*N2-N1*N2^2; % function
    FD = 3.0*(N*N1+N*N2-N1*N2);                 % first order derivative
    FDD = 6.0*N;                                % second order derivative
    Argum = 1.0-2.0*F*FDD/FD^2;
    
        if Argum>0.0
            Argum = sqrt(Argum);
        else
            Argum = 1.0;
        end
        
    DZcomb = (-2.0)*F/(FD*(1.0+Argum));         % step-length
    Zcomb = Zcomb+DZcomb;                       % new value for Zc/omega_B
end

if abs(DZcomb)<Eps
    Argum = N1*N2/(N*N1+N*N2);
    k.EoS.Omegab = (1.0-Argum)/(Zcomb-1.0);                 % Omega_B
    k.EoS.Zc = k.EoS.Omegab*Zcomb;                          % Z_c
    k.EoS.Omegaa = k.EoS.Omegab*(1.0/N-k.EoS.Omegab)*N1*N2; % Omega_A
else
    disp(' Third order method failed to calculate zcrit ')
end
