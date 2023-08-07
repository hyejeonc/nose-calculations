% Author: Oivind Wilhelmsen, Date: 2012-11-22
% This function calculates the scaled excess Gibbs energy residual 
% for the Zfac-function, to find the solution with the lowest Gibbs free 
% energy

function Gres = CB_Zfac_Gres_a(T,P,Zfac,Temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      1x1      double
%          P        Pressure                   Pa     1x1      double
%          Zfac     Compressibility factor     -      1x1      double
%          Temp     Temporary variables        -      -        struct
%
% OUTPUT:  Gres     Temporary variables        -      -        struct  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp1


if Zfac<0 
    disp(' NB: ZFac is less than zero, which is impossible ')
end

% The molar volume [m^3/kmol]
V=Zfac*ktp1.R*T/P;
% Temporart variables [m^3/mol]
N=V-Temp.b;
N1=V-ktp1.EoS.m1*Temp.b;
N2=V-ktp1.EoS.m2*Temp.b;

% Temporary variables
LNN = N2/N1;

% Check for realistic values
if LNN<0
    disp('Error, CB_Deriv is completely out of range!')
    stop
end

LNN = log(LNN);

% mole number derivative
FFN = V/N;                                                                        
     
if FFN<0
    disp('Error, CB_Deriv is completely out of range!')
    stop
end

FFN=log(FFN);                                  % (dF/dN) 

if ktp1.EoS.m1==ktp1.EoS.m2 % Van der Waals and similar EoS
    
    FF=FFN-Temp.a/(ktp1.R*T*(N1));              % F [mol]
    
else % SRK and PR and other EoS     
    
    M = ktp1.EoS.m1-ktp1.EoS.m2;
    FF=FFN-Temp.a*LNN/(M*Temp.b*ktp1.R*T);      % F [mol]
    
end

% Scaled Excess Gibbs and Helmholtz energy
Ares_RT = (FF-log(abs(Zfac)));           % [-]
Gres=Ares_RT + (Zfac-1);                 % [-] 