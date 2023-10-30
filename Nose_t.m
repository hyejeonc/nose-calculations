
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
%Parameters_animals; % default: k.Ns=24, k.Ns_art=24, k.N_cycle=1000, k.Nt=12 
%Parameters;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation starts! %%
%%%%%%%%%%%%%%%%%%%%%%%%
%function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq)
k = Parameters_animals(-30, 2, 24, 24, 4000, 0.5, 10e-6, 1, 1, 1, -2, 1.0, 1.0, 1.0, 1);
%k.tempCond = 0; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt 5: subtropical from huntley
% 6: A1 peri  7: A2 peri  8: A3 area+peri 
% 9: A4 peri  10: A5 peri 11: A6 area+peri 
% 12: Subtropical from huntely, 37 C. 
% 20: Human 


k.initCond = 0;
k.convCond = 1;
%k.initFileName = ['SA_Ns24_24_r1_1976cM_m30','.mat'];
k.fileName = ['SA_r22_4000c_m30_f10000','.mat'];
disp(k.fileName);
Nose   = CALC_OnlyAir;
save(k.fileName, 'Nose');



stop = toc
