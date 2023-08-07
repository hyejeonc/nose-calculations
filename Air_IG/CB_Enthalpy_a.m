function Enthalpy_struct = CB_Enthalpy_a(T,P,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      1x1        double
%          P        Pressure                   Pa     1x1        double
%          Z        Composition                -      1x1        double
%          phase    (1=liquid, 2=gas)          -      1x1        double
%
% OUTPUT:  Enthalpy Enthalpy                 J/kmol K    1x1     double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
global ktp1
phase   = 2;

% Z-factor, derivative and mixture parameters
[Zfac_struct, Temp]=TP_Zfac_a(T,P,Z,phase);
Zfac=Zfac_struct.Zfac;

dfdt = Temp.FFT+Temp.FFA*Temp.At;

Ures =(-ktp1.R)*T^2*(dfdt); 
Enthalpy_struct.H = Ures + ktp1.R*T*(Zfac - 1);

% Temperature derivative of the enthalpy [J/kmol K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dpdt = Temp.PT+Temp.PA*Temp.At;
d2fdt2=Temp.FFTT+Temp.FFAT*Temp.At+...
      (Temp.FFAT)*Temp.At+Temp.FFA*Temp.Att;
CVres=(-ktp1.R)*T*(T*d2fdt2+2.0*dfdt);
Enthalpy_struct.dT = CVres-T*dpdt^2/Temp.PV-ktp1.R;


% Pressure derivative of the Enthalpy [J/kmol Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Enthalpy_struct.dP = Zfac*ktp1.R*T/P + T*dpdt/Temp.PV;


% Composition derivatives of the Enthalpy [J/kmol] 
% (Partial molar enthalpy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dpdni=Temp.PN+Temp.PA*Temp.AI+Temp.PB*Temp.BI;
d2fdnidt=(Temp.FFAT)*Temp.AI+Temp.FFA*Temp.AIT+...
         (Temp.FFBT+Temp.FFAB*Temp.At).*Temp.BI;
Enthalpy_scaled=d2fdnidt+dpdni.*dpdt/(ktp1.R*T*Temp.PV)+1.0/T;
Enthalpy_struct.dN=(-Enthalpy_scaled)*ktp1.R*T^2;