% Author: Oivind Wilhelmsen, Date: 2012-11-23
% Calculates the fugacity coefficient of components in cubic EOS and 
% derivatives with respect to temperature, pressure and composition:
% ln(Phi) = mu_res/RT, F is here a part of the Helmholtz energy

function Fug_struct = TP_Fug_w(T,P,Z,phase,flag_deriv)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol       Explanation              Units    Dim.     Type
%
% INPUT:   T            Temperature                K      1 x 1    double
%          P            Pressure                   Pa     1 x 1    double
%          Z            Mole fractions             -        -      struct
%          Phase        Phase flag, 1=liq, 2=gas   -      1 x 1      -%          
%          flag_deriv   Derivative flag:           -        -      integer 
%                       0=no derivatives
%                       1=all derivatives
%
% OUTPUT:  Fug          Fugacity coeffic struct    -       -       struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
global ktp2

% Calculate the compressibility factor [-]
[Zfac_struct, Temp]=CB_Zfac_w(T,P,Z,phase,flag_deriv);
Zfac=Zfac_struct.Zfac;

% The scaled chemical potential residual [-]
DFDNI = Temp.FFN+Temp.FFA*Temp.AI+Temp.FFB*Temp.BI;
TERM = DFDNI-log(Zfac);

if TERM>30
    TERM=30*ones(size(Z));
end

% Check the calculation of the scaled Gibbs energy [-]
%X1=CB_Zfac_Gres(T,P,Zfac,Temp)
%X2=sum(TERM.*Z)

if abs(sum(Z)-1)>1E-7;
    disp(' The mole fractions do NOT sum to one, please normalize')
    sum(Z)
    stop
end    

% The fugacity coefficients Phi
Fug_struct.Fug=exp(TERM);
Fug_struct.Ln_Fug=TERM;

if flag_deriv>0
    % The temperature derivative of ln(Phi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dPdT=Temp.PT+Temp.PA*Temp.At;
    dPdNi=Temp.PN+Temp.PA*Temp.AI+Temp.PB*Temp.BI;
    d2FdNidT=Temp.FFAT*Temp.AI+Temp.FFA*Temp.AIT+...
            (Temp.FFBT+Temp.FFAB*Temp.At)*Temp.BI;
    Fug_struct.dT=d2FdNidT+dPdNi*dPdT/(ktp2.R*T*Temp.PV)+1/T;

    % The pressure derivative of ln(Phi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fug_struct.dP=-dPdNi/(ktp2.R*T*Temp.PV)-1.0/P;

    % The composition derivative of ln(Phi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IN=ones(ktp2.nc,ktp2.nc); % Large matrix with only ones

    % Diagonal matrices Ai and Bi
    BIi=(diag(Temp.BI)*IN);
    BIj=BIi';
    AIi=(diag(Temp.AI)*IN);
    AIj=AIi';
    dpdni=(diag(dPdNi)*IN);
    dpdnj=dpdni';
    
    % More caluclations 
    d2fdnidnj=(Temp.FFNB.*BIj)+(Temp.FFAB.*BIj).*AIi+Temp.FFA.*Temp.AIJ+...
              (Temp.FFNB+Temp.FFAB*AIj+Temp.FFBB*BIj).*BIi; 
               
    Fug_struct.dN=d2fdnidnj+(dpdni.*dpdnj)./(ktp2.R*T*Temp.PV)+1;
end