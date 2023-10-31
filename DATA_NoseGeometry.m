%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% extract data for nose geometry 
% input = using all data from parameter_animals
% output = Sol in "nose_geometry.mat" 
% Sol.SS = SS;
% Sol.SA = SA;
% Sol.SS_shunt = SS_shunt;
% Sol.SA_shunt = SA_shunt;
% 
% r1 : save("nose_geometry_r1.mat", "Sol");
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%% +10C Seal subtropical wo shunt, subtropical condition 
%k = Parameters_animals(10, 1, 24, 24, 1000); % Parameters_animals(tempCond, noseGeom)
%k.tempCond = 10; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
%k.fileName = ['SS_Ns24_24_r0_1000c_p10','.mat'];
%disp(k.fileName);
%Nose      = CALC_OnlyAir;
%save(k.fileName, 'Nose');




%%  0C Seal subtropical wo shunt, subtropical condition 
%k = Parameters_animals(0, 1, 24, 24, 1000); % Parameters_animals(tempCond, noseGeom)
%k.tempCond = 0; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
%k.fileName = ['SS_Ns24_24_r0_1000c_p00','.mat'];
%disp(k.fileName);
%Nose   = CALC_OnlyAir;
%save(k.fileName, 'Nose');


%% +10C Seal subtropical wo shunt, subtropical condition 
%k = Parameters_animals(10, 1, 36, 36, 1000); % Parameters_animals(tempCond, noseGeom)
%k.tempCond = 10; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
%k.fileName = ['SS_Ns36_36_r0_1000c_p10','.mat'];
%disp(k.fileName);
%disp("k.Ns")

%Nose      = CALC_OnlyAir;
%save(k.fileName, 'Nose');



%% -10C Seal subtropical wo shunt, subtropical condition 
%k = Parameters_animals(-10, 1, 24, 24, 1000); % Parameters_animals(tempCond, noseGeom)
%k.tempCond = 0; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
%k.fileName = ['SS_Ns24_24_r0_1000c_m10','.mat'];
%disp(k.fileName);
%Nose   = CALC_OnlyAir;
%save(k.fileName, 'Nose');



%% -30C Seal arctic wo shunt, arctic condition (ref)
%k = Parameters_animals(-30, 2, 24, 24, 1000); % Parameters_animals(tempCond, noseGeom)
%k.tempCond = 0; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
%k.fileName = ['SA_Ns24_24_r0_1000c_m30','.mat'];
%disp(k.fileName);
%Nose   = CALC_OnlyAir;
%save(k.fileName, 'Nose');


% -30C Seal subtropical wo shunt, arctic condition 
k1 = Parameters_animals(-30, 1, 24, 24, 1000, 0.5); % Parameters_animals(tempCond, noseGeom)
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
SS.L         = k1.L;
SS.area_a    = pchip(k1.x_a,  k1.Area_a, k1.zz);
SS.gamma_a   = pchip(k1.x_a,  k1.Per_a,  k1.zz);
SS.gamma_art = pchip(k1.x_art,k1.Per_art,k1.zz);
SS.gamma_ven = pchip(k1.x_ven,k1.Per_ven,k1.zz);


k2 = Parameters_animals(-30, 2, 24, 24, 1000, 0.5); % Parameters_animals(tempCond, noseGeom)
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
SA.L         = k2.L;
SA.area_a    = pchip(k2.x_a,  k2.Area_a, k2.zz);
SA.gamma_a   = pchip(k2.x_a,  k2.Per_a,  k2.zz);
SA.gamma_art = pchip(k2.x_art,k2.Per_art,k2.zz);
SA.gamma_ven = pchip(k2.x_ven,k2.Per_ven,k2.zz);


k3 = Parameters_animals(-30, 3, 24, 24, 1000, 0.5); % Parameters_animals(tempCond, noseGeom)
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
SS_shunt.L         = k3.L;
SS_shunt.area_a    = pchip(k3.x_a,  k3.Area_a, k3.zz);
SS_shunt.gamma_a   = pchip(k3.x_a,  k3.Per_a,  k3.zz);
SS_shunt.gamma_art = pchip(k3.x_art,k3.Per_art,k3.zz);
SS_shunt.gamma_ven = pchip(k3.x_ven,k3.Per_ven,k3.zz);


k4 = Parameters_animals(-30, 4, 24, 24, 1000, 0.5); % Parameters_animals(tempCond, noseGeom)
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
SA_shunt.L         = k4.L;
SA_shunt.area_a    = pchip(k4.x_a,  k4.Area_a,  k4.zz);
SA_shunt.gamma_a   = pchip(k4.x_a,  k4.Per_a,  k4.zz);
SA_shunt.gamma_art = pchip(k4.x_art,k4.Per_art,k4.zz);
SA_shunt.gamma_ven = pchip(k4.x_ven,k4.Per_ven,k4.zz);


Sol.SS = SS;
Sol.SA = SA;
Sol.SS_shunt = SS_shunt;
Sol.SA_shunt = SA_shunt;

save("nose_geometry_r1.mat", "Sol");
