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


% +05C Seal subtropical wo shunt 
k = Parameters_animals(-30, 2, 24, 24, 4000, 0.5, 10e-6, 1, 1, 1, 1); % Parameters_animals(tempCond, noseGeom)
%k.tempCond = 0; 
%k.noseGeom = 1;  % 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt 5: subtropical from huntley
% 6: A1 peri  7: A2 peri  8: A3 area+peri 
% 9: A4 peri  10: A5 peri 11: A6 area+peri 
% 12: Subtropical from huntely, 37 C. 

k.initCond = 0;
k.convCond = 1;

SolTemp = load('SA_r7_2000c_m30.mat');
Sol = SolTemp.Nose;

j = 18;

global k ktp1 ktp2

k = ScalesParameters(k);

gamma_a  = pchip(k.x_a,k.Per_a,k.zz);

A_a      = pchip(k.x_a,k.Area_a,k.zz);

D_a      = 4.*A_a./gamma_a;

disp(Sol.Ta(6,12))
    
for ii= 1:24
    
    xa          = [Sol.x_a(end-j,ii); (1-Sol.x_a(end-j,ii))*0.79; (1-Sol.x_a(end-j,ii))*0.21];
    disp(D_a(ii))
    disp(Sol.va(end-j,ii))
    hLam        = CALC_R_a_m_laminar(Sol.Ta(end-j,ii),xa,Sol.Tw(end-j,ii),D_a(ii),Sol.va(end-j,ii));
    hTur        = CALC_R_a_m_turbulent(Sol.Ta(end-j,ii),xa,Sol.Tw(end-j,ii),D_a(ii),Sol.va(end-j,ii));

    plotTurbulent(ii) = hTur;
    plotLaminar(ii)   = hLam;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Local entropy along z axis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nFig = nFig + 1;
figure(10); 
sz = 1.5;




%plot(zArc./L, rdNs12.plotTotSigmaQuar, '--k',   'LineWidth',sz, 'DisplayName',['N_s = 12, N_{cyc} = ', num2str(rdNs12.Ncyc/4)] ) ;hold on;
%plot(zArc./L, rdNs12.plotTotSigmaHalf, '--k',   'LineWidth',sz, 'DisplayName',['N_s = 12, N_{cyc} = ', num2str(rdNs12.Ncyc/2)] ) ;hold on;
plot(linspace(0, 0.44, 24)./0.44, plotTurbulent, ':c',   'LineWidth',sz, 'DisplayName',['Turbulent'] ) ;hold on;
plot(linspace(0, 0.44, 24)./0.44, plotLaminar,   ':k',   'LineWidth',sz, 'DisplayName',['Laminar'] ) ;hold on;

%title('Friction factor x2')
xlabel('Scaled position [-]');
ylabel('Heat transfer coefficient [W/m²K]');
set(gca,'FontSize',16);
%ylim([0 1])
set(gcf,'Position', [0 0 1080 1080]);
legend('location', 'northeast')
%saveas(figure(nFig), [pwd '/fig/reindeer/gridtest_entropy.pdf'])
%saveas(figure(nFig),'./fig/reindeer/gridtest_entropy.eps', 'epsc')
