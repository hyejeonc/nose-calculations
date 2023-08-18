function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle)
global ktp1

%%% PARAMETERS FOR CALCULATIONS FROM EXPERIMENTAL DATA
% Ambient conditions ( pressure, temperature and density)
% Relative humidity is 90 %, so multiplying 0.9 here. 
% https://www.engineeringtoolbox.com/moist-air-properties-d_1256.html 
% https://geo.libretexts.org/Bookshelves/Oceanography/Oceanography_101_(Miracosta)/08%3A_Atmospheric_Circulation/8.02%3A_Water_Moisture_in_the_Air#:~:text=The%20Water%20Cycle-,Water%20Moisture%20in%20the%20Air%20(Humidity),refers%20to%20dry%20air%20conditions.

k.pa  = 1e5;             % [Pa]
 
if tempCond == -30
    % Case 1 
    k.Ta  = 273.15-30;                % [K]
    k.xa      = 34.2/k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == -20
    % Case 2 
    k.Ta  = 273.15-20;                % [K]
    k.xa      = 92.9./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == -10 
    % Case 3 
    k.Ta  = 273.15-10;                % [K]
    k.xa      = 233.3./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == 0 
    % Case 4 
    k.Ta  = 273.15-0;                % [K]
    k.xa      = 549.7./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == 10 
    % Case 5 
    k.Ta  = 273.15+10;                % [K]
    k.xa      = 1104.3./k.pa; % 1227 * 0.9 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    
    
elseif tempCond == 20 
    % Case 5 
    k.Ta  = 273.15+20;                % [K]
    k.xa      = 2103.3./k.pa; % 2337 * 0.9 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    
end      


% Deep reindeer's body temperature
if noseGeom == 0
    k.T_body  = 273.15 + 36.3;%38.4; %***
else
    k.T_body  = 273.15 + 38.4; %***
end

k.wb      = 0.038;
k.M_dry     = ktp1.Mw(2)*0.79 + ktp1.Mw(3)*0.21;
M_vap     = ktp1.Mw(1);
k.xb      = k.wb.*k.M_dry./(M_vap.*(1-k.wb)+k.wb.*k.M_dry);

[Zfac_struct,~] = TP_Zfac_a(k.Ta,k.pa,[k.xa;0.79*(1-k.xa);0.21*(1-k.xa)],2);
rho       = Zfac_struct.c.*Ma;
rho_vap   = rho*k.wa;
k.rhoa    = rho - rho_vap;
%% DATA ABOUT REINDEER
% Reindeer mass
M_reindeer= 54;

%% ANIMAL GEOMETRY DATA 
% 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt
if noseGeom == 0 
    M_seal= 54;           %*** [kg]
    % Respiratory minute volume of air
    V_min     =  0.1; % [l/(min kg_reindeer)]
    k.V         = V_min *1e-3*(M_seal); % [m3 in 1 min]  %***
    M_tot     = k.rhoa*k.V;   % kg in 1 minuto 

else 
    M_seal= 96;%54;           %*** [kg]
    % Respiratory minute volume of air
    V_min     = 0.3967 ; %*** [l/(min kg_seal^0.75)]     0.1; % [l/(min kg_reindeer)]
    k.V         = V_min *1e-3*(M_seal^0.75); % [m3 in 1 min]  %***
    M_tot     = k.rhoa*k.V;   % kg in 1 minuto 
end


% Air mass flow
k.F_a     = M_tot/60;       % [kg/s]

disp('k.F_a === ')
disp(k.F_a)



% Blood mass flow in veins and arteries (PROPER VALUE HERE: in the whole head is 48e-4!!!!!)
k.F_b   = 2e-04*M_seal/M_reindeer;          % [kg/s] % *** A fraction of the value above, gives velocities that makes sense, see article
% Blood density
k.rho_b = 1000;             % [kg/m3]
k.rho_w = 1000;             % [kg/m3]

% Heat capacities of air, water, and blood (we assume them constant)
k.c_a   = 1003;             % [J/kg K]
k.c_w   = 4186;             % [J/kg K]
k.c_b   = 4509;            % [J/kg K] (we use the water properties for the blood)
k.k_b   = 0.505;              % We use blood conductivity both for meat and blood
k.k_w   = 0.58;
k.k_a   = 0.024;
 

%% GEOMETRICAL DATA
% 0: reindeer, 1: subtropical, 2: arctic, 3: subtropical+shunt, 4: arctic+shunt

% Lenght of the nose
k.L_reindeer       = 0.2;            % [m] %***
% Nasal cavity (area and perimeter of the cross sections) [m2 and m]
k.Area_a_reindeer  = [1500 2300 2800 2900 3200 3400 3350 3000 2950 2900 2850 2950 4150 6000 6650 5200 5350 4800 4750 4600 4500]./10.*1e-6; %*** % [m2]
k.Per_a_reindeer   = [1300 1400 1600 1800 2550 3250 4350 4000 5300 5700 5800 5500 6750 7450 5500 4500 3200 2400 1600 1100 950 ]./10.*1e-3; %***
seal = DATA_Csvdata(36, 0.5);

if noseGeom == 0 
    k.L      = k.L_reindeer;
    k.Area_a = k.Area_a_reindeer;
    k.Per_a  = k.Per_a_reindeer;
    
elseif noseGeom == 1 % 1: subtropical
    k.L      = 1e-3*(max(seal.xSubPosi) - min(seal.xSubPosi));
    k.Area_a = 1e-6*seal.xSubAreaWoshunt;
    k.Per_a  = 1e-3*seal.xSubPeri;

elseif noseGeom == 2 % 2: arctic
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*seal.xArcAreaWoshunt;
    k.Per_a  = 1e-3*seal.xArcPeri;
    
elseif noseGeom == 3 % 3: subtropical+shunt
    k.L      = 1e-3*(max(seal.xSubPosi) - min(seal.xSubPosi));
    k.Area_a = 1e-6*seal.xSubAreaWshunt;
    k.Per_a  = 1e-3*seal.xSubPeri;
    
elseif noseGeom == 4 % 4: arctic+shunt
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*seal.xArcAreaWshunt; % was Woshunt 
    k.Per_a  = 1e-3*seal.xArcPeri;
end


% result.xArc = xiArc;
% result.xArcPosi = xArcPosi; 
% result.xArcPeri = xArcPeri ;
% result.xArcAreaWshunt = xArcAreaWshunt;
% result.xArcAreaWoshunt = xArcAreaWoshunt;
% 
% result.xSub = xiSub; 
% result.xSubPosi = xSubPosi; 
% result.xSubPeri = xSubPeri ;
% result.xSubAreaWshunt = xSubAreaWshunt;
% result.xSubAreaWoshunt = xSubAreaWoshunt;

k.x_a       = linspace(0,k.L,length(k.Area_a));
k.A_a_mean  = trapz(k.x_a,k.Area_a)/k.L;
k.P_a_mean  = trapz(k.x_a,k.Per_a)/k.L;

% Mucus thickness  ?????!!!! We should find an experimental value
k.x_w     = linspace(0,k.L,length(k.Area_a));
k.d       = 0.7e-03;          % [m]
k.Per_w   = k.Per_a;
k.Area_w  = k.d.*k.Per_w;
k.A_w_mean  = trapz(k.x_w,k.Area_w)/k.L;
k.P_w_mean  = trapz(k.x_w,k.Per_w)/k.L;
k.d_mucus = 5e-4;

% Arteries (area and perimeter of the cross sections) [m2 and m]
k.x_art   = linspace(0,k.L,5);
P_art     = pchip(k.x_a,k.Per_a,k.x_art);
A_art     = [5.79 3.88 1.17 1.53 3.85]*1e-5 ;  %***max(k.Area_a)/max(k.Area_a_reindeer)    % area per m of lining 
k.Area_art= A_art.*P_art;
k.Per_art = [0.772 0.62 0.323 0.33 0.44].*P_art;
k.frac_art= [5.3 3.5 1.9 1.5 4.8]*1e-2;
k.A_art_mean  = trapz(k.x_art,k.Area_art)/k.L;
k.P_art_mean  = trapz(k.x_art,k.Per_art)/k.L;

% Veins (area and perimeter of the cross sections) [m2 and m]
k.x_ven   = linspace(0,k.L,5);
A_ven     = [6.70 6.08 1.91 2.26 5.50]*1e-4; %***       % area per m of lining 
P_ven     = pchip(k.x_a,k.Per_a,k.x_ven);
k.Area_ven= A_ven.*P_ven;
k.Per_ven = [3.31 2.76 1.70 1.51 2.20].*P_ven;
k.frac_ven= [42 38 25 20 26]*1e-2;
k.A_ven_mean  = trapz(k.x_ven,k.Area_ven)/k.L;
k.P_ven_mean  = trapz(k.x_ven,k.Per_ven)/k.L;


k.Area_tot= k.Area_art./k.frac_art;

% Meat 
k.x_m     = linspace(0,k.L,5);
k.dist_art= [0.2 0.11 0.1 0.18 0.64]*1e-03*max(k.Area_a)/max(k.Area_a_reindeer); %*** 
k.dist_ven= [0.4 0.25 0.15 0.15 0.225]*1e-03*max(k.Area_a)/max(k.Area_a_reindeer); %***
k.frac_m  = 1-k.frac_art-k.frac_ven;
k.Area_m  = k.Area_tot.*k.frac_m;
k.Area_m_mean  = trapz(k.x_m,k.Area_m)/k.L;

%% BREATHING CYCLE
if noseGeom == 0 
    k.fa      = 10;
else  
    k.fa      = 10.6;   % breaths per min %***
end

M_breath = M_tot/k.fa;    % kg/breath
k.T_breathing  = 60/k.fa/2;  % [s] (10 breaths per min)
k.A_breathing  = M_breath/2*pi/k.T_breathing;

%% DISCRETIZATION
k.Ns    = Ns;
k.Ns_art= Ns_art;
k.Nt    = 12;
k.N_cycle = N_cycle;
% Spatial discretization
k.zz    = 0:k.L/(k.Ns-1):k.L;
k.zz_art= linspace(0,k.L,k.Ns_art);

