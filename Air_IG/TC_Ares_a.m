% Author: Oivind Wilhelmsen, Date: 2012-12-19
% Calculates the residual Helmholtz energy in cubic EOS and 
% derivatives with respect to temperature, pressure and composition:
% ln(Phi) = mu_res/RT, F is here a part of the Helmholtz energy

function [Temp F_struct] =TC_Ares(T,c,Z)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol       Explanation              Units    Dim.     Type
%
% INPUT:   T            Temperature                K       1 x 1   double
%          c            Molar concentration        mol/m^3 1 x 1   double
%          Z            Composition                -       1 x nc  double
%
% OUTPUT:  Fug_stuct    Fugacity coeffic struct    -       -       struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
global ktp1

% Derivatives and parameters
Temp = CB_Mix_a(T,Z);
Temp = CB_Deriv_c_a(T,c,Temp);

% The scaled chemical potential residual [-]
dFdNi = Temp.FFN+Temp.FFA*Temp.AI+Temp.FFB*Temp.BI;
TERM = dFdNi;

if TERM>30
    TERM=30;    
end

% The scaled chemical potential residual
F_struct.mures=TERM;

% The temperature derivative of mu_res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2FdNidT=Temp.FFAT*Temp.AI+Temp.FFA*Temp.AIT+...
        (Temp.FFBT+Temp.FFAB*Temp.At)*Temp.BI;
F_struct.dT=d2FdNidT;

% The concentration derivative of mu_res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2FdNidV = Temp.FFNV+Temp.FFVA*Temp.AI+Temp.FFVB*Temp.BI;
d2FdNidc = (-1/(c^2))*d2FdNidV;

F_struct.dc=d2FdNidc;

% The volume derivative:
F_struct.dV=d2FdNidV;

% The composition derivative of mu_res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IN=ones(ktp1.nc,ktp1.nc); % Large matrix with only ones

% Diagonal matrices Ai and Bi
BIi=(diag(Temp.BI)*IN);
BIj=BIi';
AIi=(diag(Temp.AI)*IN);
AIj=AIi';
d2FdNidV_mat=(diag(d2FdNidV)*IN);  % The d2F_dNidV matrix, since (V/n)
                                    % depends on ni.
                                    
% More caluclations 
d2fdnidnj=(Temp.FFNB.*BIj)+(Temp.FFAB.*BIj).*AIi+Temp.FFA.*Temp.AIJ+...
          (Temp.FFNB+Temp.FFAB*AIj+Temp.FFBB*BIj).*BIi;  % dln(Phi)/dNi
        
               
Matrix=d2fdnidnj+(1/c)*d2FdNidV_mat;     % dln(Phi)/dV*dV/dNi;
F_struct.dN=Matrix;                      % Constant C
F_struct.dN_v=d2fdnidnj;                 % Constant V
