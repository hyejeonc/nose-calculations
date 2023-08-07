% Author: Oivind Wilhelmsen, Date: 2012-11-22
% Calculates the mixture parameters in cubic EoS and the derivatives
% with respect to composition and temperature 
% 2010-01-13, compared to TPLib, gives equal results

function Temp = CB_Mix_w(T,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K     1 x 1     double
%          Z        Mole fractions             -     ncx 1     double
%
% OUTPUT:  Temp     Temporary variables, such
%                   as derivatives and mixing 
%                   parameters, A and B        -      (-)      struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ktp2 

% Calculate Fa and Fb
Fa=ktp2.EoS.Omegaa.*(ktp2.R*ktp2.Tc).^2./ktp2.Pc;
Fb=ktp2.EoS.Omegab.*ktp2.R*ktp2.Tc./ktp2.Pc;

% More parameters
S=ktp2.EoS.alfa+ktp2.EoS.beta*ktp2.Acf-ktp2.EoS.gamma.*ktp2.Acf.^2;
XX=S.^2;
YY=-2.0.*S.*(S+1.0);

% Reduced temperatures [-]
Tr=T./ktp2.Tc;

% Alpha values [-]
alpha=1+XX.*(Tr-1)+YY.*(sqrt(Tr)-1);

% Alpha combined value [-]
afac=alpha.*Fa;

% Variable
sqrt_afac2=sqrt(afac*afac');

% Alpha matrix with interaction parameters [-]    
alpha_matrix=sqrt_afac2.*(1-ktp2.EoS.kij);

% The parameter a [Jm^3/kmol^2]
Temp.a=sum(sum((Z*Z').*alpha_matrix));

% The parameter b [m^3/kmol]
Temp.b=sum(Z.*Fb);


% (1/N)dA_s/dNi
Temp.AI=2*alpha_matrix*Z;

% dai_dT
dai_dT=Fa.*((XX./ktp2.Tc)+0.5*(YY)./sqrt(T*ktp2.Tc));

% d2ai_dT2
d2ai_dT2=-Fa.*(YY)./(4*sqrt(ktp2.Tc)*T^(1.5));

% daij_dT
Inb2=afac*dai_dT';
daij_dT=(1-ktp2.EoS.kij).*(1./(2.*sqrt_afac2)).*(Inb2+(Inb2'));

% In between variable:
Inb=((afac*dai_dT'+(afac*dai_dT')').^2)./(sqrt_afac2.^2);

% d2aij_dT2
d2aij_dT2=0.5*(1-ktp2.EoS.kij).*(1./sqrt_afac2).*(2*dai_dT*dai_dT'+afac*d2ai_dT2'+(afac*d2ai_dT2')'-0.5*(Inb));

% (1/N)d^2A_s/dNidT
Temp.AIT=2*daij_dT*Z;

% d2A/dNidNj
Temp.AIJ=2*alpha_matrix;

% (1/N^2)dA/dT
Temp.At=0.5*sum(Z'*Temp.AIT);

% (1/N^2)dA2/dT2
Temp.Att=sum(sum(((Z*Z').*d2aij_dT2)));

% dB/dNi
Temp.BI=Fb;