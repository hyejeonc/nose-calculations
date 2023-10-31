%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% main running file for seal nose. 
% r0 : default revision. medium direction in vein is -1 during exhalation (dss020) 
% r1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

global k ktp1 ktp2
tic;
% Add folders to the path
WorkDirectory=pwd;                  
cd ..;cd Water_EoS;WTBasisDirectory=pwd;cd(WorkDirectory);
cd ..;cd Air_IG;ATBasisDirectory=pwd;cd(WorkDirectory);

addpath(WTBasisDirectory)   % Thermo basis
addpath(ATBasisDirectory)   % Thermo basis

% Select the components:
% C1 C2 C3 iC4 nC4 iC5 nC5 nC6 nC7 nC8 nC9 H2O N2 CO2 H2S cC6 Ar O2 H2
% 1  2  3  4   5   6   7   8   9   10  11  12  13 14  15  16  17 18 19

Comp1=[12 13 18];
Comp2=[12];

% Select the Equation of State
k.EoS = 3;                  % 1 = Van der Waals EoS
                            % 2 = Soave Redlich Kwong EoS
                            % 3 = Peng Robinson EoS

% Parameter for the thermodinamic package
ktp1 = TP_Init(Comp1,k.EoS);
ktp2 = TP_Init(Comp2,k.EoS); 


Ta = 240 ;
pa = 100000 ;
xa = [0.0595; 0.7430; 0.1975]; 
Mw = [18.015, 28.014, 44.010];
Ma = sum(xa.*Mw);

Tb = 310; 

[Zfac_struct,~] = TP_Zfac_a(Ta, pa , xa ,2);
disp("** == Ta = 240 == ** ")
disp("== air == ")
disp("== Zfac.c == ")
disp(Zfac_struct.c)

mu_a        = TP_ChemicalPotential_a(Ta, pa ,xa, 1);
disp("== mu_a == ")
disp(mu_a)

h_id        = IG_H_a(Ta, xa);
h_res       = CB_Enthalpy_a(Ta, pa, xa);
h_i         = (h_id.dZ + h_res.dN)./Mw;
h_a         = (h_id.H + h_res.H)./Ma;
disp("== h_id == ")
disp(h_id)
disp("== h_res == ")
disp(h_res)
disp("== h_i == ")
disp(h_i)
disp("== h_a == ")
disp(h_a)


disp("== water == ")
mu_m        = TP_ChemicalPotential_w(Ta,pa,1);
disp("== mu_m == ")
disp(mu_m)

h_wa = h_i(1);
h_id        = IG_H_w(Ta); % water enthalpy 
h_res       = CB_Enthalpy_w(Ta, pa);
h_w         = (h_id.H + h_res.H)./Mw;
%X           = [1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii); -(mu_m.mu(1)./Sol.Tw(j,ii)-mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii))];

disp("== h_wa == ")
disp(h_wa)
disp("== h_id == ")
disp(h_id)
disp("== h_res == ")
disp(h_res)
disp("== h_w == ")
disp(h_w)s



disp("** == Tb = 310 == ** ")
Ta = Tb;

[Zfac_struct,~] = TP_Zfac_a(Ta, pa , xa ,2);
%disp("** == Ta = 240 == ** ")
disp("== air == ")
disp("== Zfac.c == ")
disp(Zfac_struct.c)

mu_a        = TP_ChemicalPotential_a(Ta, pa ,xa, 1);
disp("== mu_a == ")
disp(mu_a)

h_id        = IG_H_a(Ta, xa);
h_res       = CB_Enthalpy_a(Ta, pa, xa);
h_i         = (h_id.dZ + h_res.dN)./Mw;
h_a         = (h_id.H + h_res.H)./Ma;
disp("== h_id == ")
disp(h_id)
disp("== h_res == ")
disp(h_res)
disp("== h_i == ")
disp(h_i)
disp("== h_a == ")
disp(h_a)


disp("== water == ")
mu_m        = TP_ChemicalPotential_w(Ta,pa,1);
disp("== mu_m == ")
disp(mu_m)

h_wa = h_i(1);
h_id        = IG_H_w(Ta); % water enthalpy 
h_res       = CB_Enthalpy_w(Ta, pa);
h_w         = (h_id.H + h_res.H)./Mw;
%X           = [1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii); -(mu_m.mu(1)./Sol.Tw(j,ii)-mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii))];

disp("== h_wa == ")
disp(h_wa)
disp("== h_id == ")
disp(h_id)
disp("== h_res == ")
disp(h_res)
disp("== h_w == ")
disp(h_w)            
            
% h_tot(1,ii) = (h_id.H + h_res.H)./Ma(ii);
% h_i(:,ii)   = (h_id.dZ + h_res.dN)./ktp1.Mw;
% h_a(1,ii)   = h_i(1,ii);
% mu_str      = TP_ChemicalPotential_w(Tm(ii),k.pa,1);
% mu_m(1,ii)  = mu_str.mu;
% h_id        = IG_H_w(Tm(ii));
% h_res       = CB_Enthalpy_w(Tm(ii),k.pa);
% h_w(1,ii)   = (h_id.H + h_res.H)./ktp2.Mw;            

disp("** == Ta = 240 == ** ")
#disp("** == xa = [0.1, 0.7, 0.2] == ** ")
xa = [0.3; 0.5; 0.2]; 
Ta = 240 ;

[Zfac_struct,~] = TP_Zfac_a(Ta, pa , xa ,2);
%disp("** == Ta = 240 == ** ")
disp("== air == ")
disp("== Zfac.c == ")
disp(Zfac_struct.c)

mu_a        = TP_ChemicalPotential_a(Ta, pa ,xa, 1);
disp("== mu_a == ")
disp(mu_a)

h_id        = IG_H_a(Ta, xa);
h_res       = CB_Enthalpy_a(Ta, pa, xa);
h_i         = (h_id.dZ + h_res.dN)./Mw;
h_a         = (h_id.H + h_res.H)./Ma;
disp("== h_id == ")
disp(h_id)
disp("== h_res == ")
disp(h_res)
disp("== h_i == ")
disp(h_i)
disp("== h_a == ")
disp(h_a)


disp("== water == ")
mu_m        = TP_ChemicalPotential_w(Ta,pa,1);
disp("== mu_m == ")
disp(mu_m)

h_wa = h_i(1);
h_id        = IG_H_w(Ta); % water enthalpy 
h_res       = CB_Enthalpy_w(Ta, pa);
h_w         = (h_id.H + h_res.H)./Mw;
%X           = [1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii); -(mu_m.mu(1)./Sol.Tw(j,ii)-mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii))];

disp("== h_wa == ")
disp(h_wa)
disp("== h_id == ")
disp(h_id)
disp("== h_res == ")
disp(h_res)
disp("== h_w == ")
disp(h_w)         

