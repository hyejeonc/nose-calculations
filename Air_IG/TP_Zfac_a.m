function [Zfac_struct, Temp]=TP_Zfac_a(T,P,Z,phase)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation              Units    Dim.     Type
%
% INPUT:   T        Temperature                K      1 x 1    double
%          P        Pressure                   Pa     1 x 1    double
%          Z        Mole fractions             -    noc x 1    double
%          phase    Phase(1=liquid, 2=gas)     -      1 x 1    double
%
% OUTPUT:  Zfac     Compressibility factor     -     1x1       struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
global ktp1
    
    % Calculate the parameters from mixing-rules
    Temp = CB_Mix_a(T,Z);
    
    % Scaled parameter A [-]
    A_s=Temp.a*P/((T*ktp1.R)^2);
    
    % Scaled parameter B [-]
    B_s=Temp.b*P/(T*ktp1.R);
    
    % Parameters
    m1_m2_p=(ktp1.EoS.m1+ktp1.EoS.m2);
    m1_m2_m=(ktp1.EoS.m1*ktp1.EoS.m2);
    
    PP = -B_s*(m1_m2_p)-B_s-1;
    QQ =  (B_s^2)*m1_m2_m+(B_s^2)*(m1_m2_p)+B_s*m1_m2_p+A_s;
    RR = -1*(B_s^3)*m1_m2_m-B_s.^2*m1_m2_m-A_s*B_s;
    
    % The coefficients of the polynomial Z^3+PP*Z^2+QQ*Z+RR=0, where
    % Z is the compressibility
    p=[1 PP QQ RR];
    
    % The solution of the polynomial
    r_all=roots(p);
    
    % Extract only the real roots
    r_real=r_all(imag(r_all)==0);
    
    if length(r_real)>1  % More than one root
        
        Zfac_l=min(r_real);  % Liquid
        Zfac_g=max(r_real);  % Gas
        
        
        if phase==0    % Minimum Gibbs energy
            
            % Gibbs energy residuals [J/kmol]
            G_res_g=CB_Zfac_Gres_a(T,P,Zfac_g,Temp);
            G_res_l=CB_Zfac_Gres_a(T,P,Zfac_l,Temp);
            
            if G_res_g<=G_res_l % Gas is more stable than liquid
                Zfac=Zfac_g;
            else                % Liquid is more stable than gas
                Zfac=Zfac_l;
            end
            
        elseif phase==1
            Zfac=Zfac_l;  % Liquid
        elseif phase==2
            Zfac=Zfac_g;  % Gas
        else
            disp('The wrong phase has been chosen')
        end
    else    % Only one root
        Zfac=r_real;
    end
    
    % HERE WE SHOULD ADD SOMETHING TO INCREASE THE ACCURACY OF ZFAC
    
    % Calculate derivatives of the Helmholz energy
    Temp = CB_Deriv_a(T,P,Zfac,Temp);
    
    % Temperature derivative of compressibility
    dPdT=Temp.PT+Temp.PA*Temp.At;         % (dP/dT) [Pa/K]
    dVdT=-dPdT/Temp.PV;                   % (dV/dT)_{p,Ni} [m^3/kmol K]
    Zfac_dT=P*dVdT/(ktp1.R*T)-Zfac/T;      % (dZfac/dT)_{p,Ni}
    
    % Pressure derivative of the compressibility
    Zfac_dP=P/(Temp.PV*ktp1.R*T) + Zfac/P;   % (dZfac/dP)_{T,Ni}
    
    % Composition derivative of the compressibility
    DPDNI=Temp.PN+Temp.PA*Temp.AI+...     % N_tot*(dP/dNi)
        Temp.PB*Temp.BI;
    DVDNI=-DPDNI/Temp.PV;                 % (dV/dNi)
    Zfac_dN=P/(ktp1.R*T)*DVDNI - Zfac;     % Scaled compressibility derivative
    % (1/N_tot)*(dZfac/dNi)
    
    % Allocate values to struct
    Zfac_struct.Zfac=Zfac;
    Zfac_struct.c=P/(Zfac*ktp1.R*T);       % Concentration [mol/m^3]
    Zfac_struct.dT=Zfac_dT;
    Zfac_struct.dP=Zfac_dP;
    Zfac_struct.dN=Zfac_dN;
