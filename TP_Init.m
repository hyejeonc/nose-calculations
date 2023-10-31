function k = TP_Init(Comp,EoS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Symbol   Explanation                   Dimensions    Type
%
% INPUT:   Comp     Index of components             1 x nc.     double
%          EoS      Equation of state:              1 x 1       double
%                   0: Van der Waals 
%                   1: Soave Redlich Kwong
%                   2: Peng Robinson
%
% OUTPUT:  k        Struct containing constants       -        struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Components
k.Comp = Comp;

% Tolerances:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rachford Rice tolerances
k.Tol_RR1=1e-6;
k.Tol_RR2=1e-6;

% Change in fugacities in nested loop;
k.Tol_XY=1e-7;
k.MAX_iter_RR=400;
k.MAX_iter_RR_NR=30;

% Maximum number of iterations for saturation-computations
k.MAX_iter_Sat=120;
k.MAX_iter_Sat_NEWTON=30;

% Tolerance in stability equations
k.TOL_SS=1e-4;

% Other choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution flags
k.Flag_print=false;      % Print flag
k.Flag_RR_DEM=true;     % Use DEM in RR-calculation (Every 5th iteration)
k.Flag_RR_NR=true;      % Use Second order approach (After 11 iterations)    



% Components:
% C1 C2 C3 iC4 nC4 iC5 nC5 nC6 nC7 nC8 nC9 H2O N2 CO2 H2S c-C6 Ar O2 H2
% 1  2  3  4   5   6   7   8   9   10  11  12  13 14  15  16   17 18 19

% Critical temperatures [K]
Tc=[190.56
    305.32
    369.83
    407.85
    425.12
    460.39
    469.70
    507.60
    540.20
    568.70
    594.60
    647.14
    126.20
    304.12
    373.40
    553.5
    150.86
    154.58
    32.98];

% Critical pressures [Pa]
Pc=[45.99
    48.72
    42.48
    36.40
    37.96
    33.81
    33.70
    30.25
    27.40
    24.90
    22.90
    220.64
    33.98
    73.74
    89.63
    40.73
    48.98
    50.43
    12.93]*1E5;

% Pitzers acentric factor [-]
Acf=[0.011
    0.099
    0.152
    0.186
    0.200
    0.229
    0.252
    0.300
    0.350
    0.399
    0.445
    0.344
    0.037
    0.225
    0.090
    0.211
   -0.002
    0.021
   -0.217];

% The reference enthalpy [J/kmol]
Href=[-7.452
    -8.382
    -10.468
    -13.499
    -12.579
    -15.37
    -14.676
    -16.694
    -18.765
    -20.875
    -22.874
    -24.1814
    0.00
    -39.351
    -2.063
    -12.33
    0.00
    0.00
    0.00]*1E7;

% The reference entropy [J/K kmol]
Sref=[1.8627
    2.2912
    2.702
    2.955
    3.0991
    3.4374
    3.4945
    3.8874
    4.2798
    4.6723
    5.064
    1.88724
    1.915
    2.13677
    2.056
    2.97276
    1.54737
    2.05043
    1.30571]*1E5;

% Coefficients of ideal gas heat capacity polynomial: 
Cpi=[4.568	-8.975	3.631	-3.407	1.091
    4.178	-4.427	5.660	-6.651	2.487
    3.847	5.131	6.011	-7.893	3.079
    3.351	17.883	5.477	-8.099	3.243
    5.547	5.536	8.057	-10.571	4.134
    1.959	38.191	2.434	-5.175	2.165
    7.554	-0.368	11.846	-14.939	5.753
    8.831	-0.166	14.302	-18.314	7.124
    9.634	4.156	15.494	-20.066	7.770
    10.824	4.983	17.751	-23.137	8.980
    12.152	4.575	20.416	-26.770	10.465
    4.395	-4.186	1.405	-1.564	0.632
    3.539	-0.261	0.007	0.157	-0.099
    3.259	1.356	1.502	-2.374	1.056
    4.266	-3.438	1.319	-1.331	0.488
    4.035	-4.433	16.834	-20.775	7.746
    2.500	0.000	0.000	0.000	0.000
    3.630	-1.794	0.658	-0.601	0.179
    2.883	3.681	-0.772	0.692	-0.213];

% Translate Cp-coefficient values
Cp(:,1)=Cpi(:,1);
Cp(:,2)=Cpi(:,2)*1E-3;
Cp(:,3)=Cpi(:,3)*1E-5;
Cp(:,4)=Cpi(:,4)*1E-8;
Cp(:,5)=Cpi(:,5)*1E-11;

% Minimum temperature of Cp-polynomial [K]
Cp_Tmin=[50
    50
    50
    50
    200
    200
    200
    200
    200
    200
    200
    50
    50
    50
    50
    100
    50
    50
    50];

% Maximum temperature of Cp-polynomial [K]
Cp_Tmax=[1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000
    1000];

% The moleweights [kg/kmol]
Mw=[16.043
    30.070
    44.097
    58.123
    58.123
    72.150
    72.150
    86.177
    100.204
    114.231
    128.258
    18.015
    28.014
    44.010
    34.082
    84.161
    39.948
    31.999
    2.016];


% Extract the necessary variables
k.nc=length(Comp);      % number of components
k.Tc=Tc(Comp);          % critical temperature [K]
k.Pc=Pc(Comp);          % critical pressure [Pa]
k.Acf=Acf(Comp);        % acentric factor [-]
k.Mw=Mw(Comp);          % moleweight [kg/kmole]
k.H_std=Href(Comp);     % Standard enthalpy [J/kmole]
k.S_std=Sref(Comp);     % Standard entropy [J/kmole K]
k.Cp_Tmin=Cp_Tmin(Comp);% Minimum temperature cp [K]
k.Cp_Tmax=Cp_Tmax(Comp);% Minimum temperature cp [K]


for i=1:k.nc   
    k.Cp(i,:)=Cp(Comp(i),:); % Heat capacity coefficients
end

% The universal gas constant [J/kmol K]
k.R=8.314462175*1E3;
% Avogadro number [1/kmol]
k.Na=6.0221415*1E26;

% Parameters needed for the Equation of State (EoS)

if    EoS==1               % Van der Waals EoS

      k.EoS.m1=0.0;
      k.EoS.m2=0.0;
      k.EoS.alfa=0.0;
      k.EoS.beta=0.0;
      k.EoS.gamma=0.0; 
      
      k.EoS.kij=zeros(k.nc,k.nc);  % Interaction parameters
      
elseif EoS==2               % Soave Redlich Kwong EoS
    
      k.EoS.m1=-1.0;
      k.EoS.m2=0.0;
      k.EoS.alfa=0.48;
      k.EoS.beta=1.574;
      k.EoS.gamma=0.176;  
      
      k.EoS.kij=zeros(k.nc,k.nc);  % Interaction parameters
                
elseif EoS==3 % Peng Robinson EoS    
    
      k.EoS.m1=-1.0-sqrt(2);
      k.EoS.m2=-1.0+sqrt(2);
      k.EoS.alfa=0.37464;
      k.EoS.beta=1.54226;
      k.EoS.gamma=0.26992;
      
      k.EoS.kij=zeros(k.nc,k.nc);  % Interaction parameters
          
else 
    disp('This EoS is unfortunately not yet available')
end


% Calculate Omega_A and Omega_B for the EoS:
k = CB_constants(k,EoS);

% The reference temperature and pressures:
k.Tref=298.15;   % [K]
k.Pref=1E5;      % [Pa]