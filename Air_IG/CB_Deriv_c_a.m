% Author: Oivind Wilhelmsen, Date: 2012-11-27
% This function calculates the derivatives necessary for the rest of the
% state functions in cubic equations of state, based on the concentration

function Temp = CB_Deriv_c(T,c,Temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K       1x1     double
%          c        Concentration              mol/m^3 1x1     double
%          Temp     Temporary variables        -       -       struct
%
% OUTPUT:  Temp     Temporary variables        -       -       struct  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1

% The molar volume [m^3/kmol]
V=1/c;

% Temporart variables [m^3/mol]
N=V-Temp.b;
N1=V-ktp1.EoS.m1*Temp.b;
N2=V-ktp1.EoS.m2*Temp.b;

% Temporary variables
DEN=N1*N2;
LNN = N2/N1;

LNN = log(LNN);
M = ktp1.EoS.m1-ktp1.EoS.m2;

% Parameters needed in subsequent derivatives

%%%% Pressure derivatives %%%%%%
Temp.PN = ktp1.R*T/N;                                   % N_tot*(dP/dN_tot)
Temp.PA = (-1.0)/DEN;                                  % N_tot^2*(dP/dA)
Temp.PB = ktp1.R*T/N^2-Temp.a*...                       % N_tot*(dP/dB)
         (ktp1.EoS.m1*N2+ktp1.EoS.m2*N1)/DEN^2;  
Temp.PT = ktp1.R/N;                                     % (dP/dT)
Temp.PV = (-ktp1.R)*T/N^2+Temp.a*(N1+N2)/DEN^2;         % N_tot*(dP/dV)

% Derivatives of the scaled Helmholtz energy resiudal FF. 
Temp.FFN = V/N;

%if Temp.FFN<0
%    disp('Error, CB_Deriv is completely out of range!')
%    save('Dump')
%    stop
%end

Temp.FFNV=-Temp.b/(V*N);                                % (d^2F/dVdN)
Temp.FFNB = 1.0/N;                                      % (d^2F/dNdB) 
Temp.FFN=log(Temp.FFN);                                 % (dF/dN) 

if ktp1.EoS.m1==ktp1.EoS.m2 % Van der Waals and similar EoS
    
    Temp.FF=Temp.FFN-Temp.a/(ktp1.R*T*(N1));              % F [mol]
    Temp.FFA=-1/(ktp1.R*T*N1);                            % N_tot*(dF/dA)
    Temp.FFB=1.0/N-Temp.a*ktp1.EoS.m1/(ktp1.R*T*N1^2);     % (dF/dB)
    Temp.FFT = Temp.a/(ktp1.R*T^2*N1);                    % (1/N_tot)*(dF/dT)
    Temp.FFAB = -ktp1.EoS.m1/(ktp1.R*N1^2);                % N_tot^2*(d^2F/dAdB)
    Temp.FFAT = 1/(ktp1.R*T^2*N1);                        % N_tot*(d^2F/dAdT)
    Temp.FFBB = 1/N^2-2*Temp.a*ktp1.EoS.m1^2/(ktp1.R*T*N1^2);  % N_tot*(d^2F/dB^2)          
    Temp.FFBT = Temp.a*ktp1.EoS.m1/(ktp1.R*T^2*N1^2);      % (d^2F/dBdT)
    Temp.FFTT=(-2.0)*Temp.a/(ktp1.R*T^3*N1);              % (1/N_tot)*(d^F/dT^2)    
    Temp.FFV=-Temp.b/(V*N)+Temp.a/(ktp1.R*T*N1^2);        % (1/N_tot)*(dF/dV)
    Temp.FFVT=-Temp.a/(ktp1.R*T^2*N1^2);                  % (1/N_tot)*(d^2
    Temp.FFVA=1/(ktp1.R*T*N1^2);                          % N_tot*(d^2F/dVdA)
    Temp.FFVV=-1/V^2+1/N^2-2*Temp.a/(ktp1.R*T*N1^2);      % (1/N_tot)*(d^2F/dV^2)
    Temp.FFVB=-1/N^2+2*Temp.a*ktp1.EoS.m1/(ktp1.R*T*N1^3); % (d^2F/dVdB)
    
else % SRK and PR and other EoS
    
    Temp.FF=Temp.FFN-Temp.a*LNN/(M*Temp.b*ktp1.R*T);      % F [mol]
    Temp.FFA=(-LNN)/(M*Temp.b*ktp1.R*T);                  % N_tot*(dF/dA)
    Temp.FFB=1.0/N-Temp.a/(Temp.b*ktp1.R*T)*V/DEN+ ...    % (dF/dB)
            Temp.a/(Temp.b^2*ktp1.R*T*M)*LNN;
    Temp.FFT = Temp.a/(M*Temp.b*ktp1.R*T^2)*LNN;          % (1/N_tot)*(dF/dT)
    Temp.FFAB = LNN/(M*Temp.b^2*ktp1.R*T)-...             % N_tot^2*(d^2F/dAdB)
                V/(Temp.b*ktp1.R*T*DEN);
    Temp.FFAT = LNN/(M*Temp.b*ktp1.R*T^2);                % N_tot*(d^2F/dAdT)
    Temp.FFBB = 1.0/N^2-2.0*Temp.a*LNN/(Temp.b^3*ktp1.R*T*M)+ ...
              Temp.a*V/(ktp1.R*T*(Temp.b*DEN)^2)*...
             (DEN-Temp.b*(ktp1.EoS.m1*N2+ ...
              ktp1.EoS.m2*N1))+...
             Temp.a*V/(Temp.b^2*ktp1.R*T*DEN);            % N_tot*(d^2F/dB^2)          
    Temp.FFBT = Temp.a*V/(Temp.b*ktp1.R*T^2*DEN)-...      % (d^2F/dBdT)
            Temp.a*LNN/(ktp1.R*(Temp.b*T)^2*M);
    Temp.FFTT=(-2.0)*Temp.a*LNN/(M*Temp.b*ktp1.R*T^3);    % (1/N_tot)*(d^F/dT^2)
    Temp.FFV=-Temp.b/(V*N)+Temp.a/(N1*N2*ktp1.R*T);       % (1/N_tot)*(dF/dV)
    Temp.FFVT=-Temp.a/(N1*N2*ktp1.R*T^2);                 % (1/N_tot)*(d^2
    Temp.FFVA=1/(N1*N2*ktp1.R*T);                         % N_tot*(d^2F/dVdA)
    Temp.FFVV=Temp.b*(2.0*V-Temp.b)/(V^2*N^2)- ...       % (1/N_tot)*(d^2F/dV^2)
        (Temp.a/(ktp1.R*T))*((N1+N2)/(N1^2*N2^2));
    Temp.FFVB=-1/N^2+(Temp.a/(ktp1.R*T))*...
        (ktp1.EoS.m2/(N1*(N2^2))+ktp1.EoS.m1/((N1^2)*N2)); % (d^2F/dVdB)
end