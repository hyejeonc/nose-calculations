% Author: Oivind Wilhelmsen, Date: 2012-11-22
% Calculates the resiudal Entropy of cubic EOS and the 
% derivatives with respect to temperature, pressure and composition
% The entropy is defined as the derivative of the helmholz energy, A
% (dA/dT)_{v,ni}=S;

function Entropy_struct =CB_Entropy_a(T,P,Z)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      1x1        double
%          P        Pressure                   Pa     1x1        double
%          Z        Composition                -      1x1        double
%          phase    (0=min Gibbs energy
%                    1=liquid, 2=gas)          -      1x1        double
%
% OUTPUT:  Entropy  Entropy                  J/kmol K    1x1     double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
global ktp1
phase = 2;
% Z-factor, derivative and mixture parameters
[Zfac_struct Temp]=TP_Zfac_a(T,P,Z,phase);
Zfac=Zfac_struct.Zfac;

dPdT = Temp.PT+Temp.PA*Temp.At;
dFdT = Temp.FFT+Temp.FFA*Temp.At;

% The entropy [J/kmole K]
Entropy_struct.S=(-ktp1.R)*(T*dFdT+Temp.FF-log(abs(Zfac)));


% The entropy temperature derivative [J/kmol K^2]
d2FdT2=Temp.FFTT+Temp.FFAT*Temp.At+(Temp.FFAT)*Temp.At+...
       Temp.FFA*Temp.Att;
CVres = (-ktp1.R)*T*(T*d2FdT2+2.0*dFdT);
Entropy_struct.dT= CVres/T-dPdT^2/Temp.PV-ktp1.R/T;


% The entropy pressure derivative [J/kmol K Pa]
Entropy_struct.dP= dPdT/Temp.PV + ktp1.R/P;


% The entropy composition derivative [J/kmol K]
dPdNi=Temp.PN+Temp.PA*Temp.AI+Temp.PB*Temp.BI;
d2FdNidT=(Temp.FFAT)*Temp.AI+Temp.FFA*Temp.AIT+...
         (Temp.FFBT+Temp.FFAB*Temp.At)*Temp.BI;

% The residual enthalpy [J/kmol]
Enthalpy_scaled=d2FdNidT+dPdNi*dPdT/(ktp1.R*T*Temp.PV)+1.0/T;
Enthalpy=(-Enthalpy_scaled)*ktp1.R*T^2;

% dHresdT is now in placed in res
dFdNi = Temp.FFN+Temp.FFA*Temp.AI+Temp.FFB*Temp.BI;  % (dA/dNi)
Term = dFdNi-log(abs(Zfac));                         % log(Fug(i))
    
Entropy_struct.dN=Enthalpy/T-ktp1.R*Term;