function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, vLimit)
global ktp1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% main running file for seal nose. 
% r0 : default revision. medium direction in vein is -1 during exhalation (dss020) 
% r1 : correct signs in DIFF_OnlyAir.m, DIFF_OnlyAir2.m, CALC_OnlyAir.m
% r2 : body temperature to 35.8 C for subtropical due to validation with Huntley 1984 (it was 37 C)
% r3 : body temp 35.8  and body temp 36.3 from Huntley and Folkow 
% r4 : change area and perimeter flipped
% r5 : change area-area 
% r6 : change mucus thickness to 10 um.  
% r7 : add resistance for R-am, increase. 
% r8 : add resistance for R-ij, decrease. variable about factor for resistance (conductivity).
% r9 : add resistance for R-am, all elements increased. 
% r10: add friction factor.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r20: human nose
% r21: seal nose w/ Vtide (tidal volume) and Freq (respiration frequency) as variables 
% r22: r21 was with c_b 4201 not 3XXX, so modify this. additionally, study about laminar flow correlation also. 


%%% PARAMETERS FOR CALCULATIONS FROM EXPERIMENTAL DATA
% Ambient conditions ( pressure, temperature and density)
% Relative humidity is 90 %, so multiplying 0.9 here. 
% https://www.engineeringtoolbox.com/moist-air-properties-d_1256.html 
% https://geo.libretexts.org/Bookshelves/Oceanography/Oceanography_101_(Miracosta)/08%3A_Atmospheric_Circulation/8.02%3A_Water_Moisture_in_the_Air#:~:text=The%20Water%20Cycle-,Water%20Moisture%20in%20the%20Air%20(Humidity),refers%20to%20dry%20air%20conditions.

k.pa  = 1e5;             % [Pa]

if tempCond == -40
    % Case 1 
    k.Ta  = 273.15-40;                % [K] 
    k.xa      = 11.6/k.pa;    % 12.84 * 0.9 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;

elseif tempCond == -35
    % Case 1 
    k.Ta  = 273.15-35;                % [K]
    k.xa      = 22.5/k.pa;  % 25 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
     
elseif tempCond == -30
    % Case 1 
    k.Ta  = 273.15-30;                % [K]
    k.xa      = 34.2/k.pa; % 38 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
    disp('k.xa === ')
    disp(k.xa)
    disp('Ma === ')
    disp(Ma)
    disp('k.wa === ')
    disp(k.wa)
     
     
elseif tempCond == -3020
	% Tamb = -30 C, RH = 20 %
    k.Ta  = 273.15-30;                % [K]
    k.xa      = 7.6/k.pa; % 38 * 0.2
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
       
elseif tempCond == -3080 % RH = 80 % 
	% Tamb = -30 C, RH = 80 %
    k.Ta  = 273.15-30;                % [K]
    k.xa      = 30.4/k.pa; % 38 * 0.8
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == -25
    % Case 2 
    k.Ta  = 273.15-25;                % [K]
    k.xa      = 56.9./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
       
elseif tempCond == -20
    % Case 2 
    k.Ta  = 273.15-20;                % [K]
    k.xa      = 92.9./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;

elseif tempCond == -15
    % Case 3 
    k.Ta  = 273.15-15;                % [K]
    k.xa      = 148.7./k.pa;  % 168.2 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == -10 
    % Case 3 
    k.Ta  = 273.15-10;                % [K]
    k.xa      = 233.3./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == -5
    % Case 3 
    k.Ta  = 273.15-5;                % [K]
    k.xa      = 361.35./k.pa;  % 401.5*0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == 0 
    % Case 4 
    k.Ta  = 273.15-0;                % [K]
    k.xa      = 549.7./k.pa;
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;

elseif tempCond == 5 
    % Case 5 **
    k.Ta  = 273.15+5;                % [K]
    k.xa      = 784.71./k.pa; % 871.9 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == 10 
    % Case 6 
    k.Ta  = 273.15+10;                % [K]
    k.xa      = 1104.3./k.pa; % 1227 * 0.9 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    
     
elseif tempCond == 1080 % RH = 80 % 
    % Case 6 
    k.Ta  = 273.15+10;                % [K]
    k.xa      = 981.6./k.pa; % 1227 * 0.8 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    

elseif tempCond == 1020 % RH = 20 % 
    % Case 6 
    k.Ta  = 273.15+10;                % [K]
    k.xa      = 245.4./k.pa; % 1227 * 0.2 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    
    
elseif tempCond == 15 
    % Case 7 **
    k.Ta  = 273.15+15;                % [K]
    k.xa      = 1533.6./k.pa; % 1704 * 0.9 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    
    
elseif tempCond == 20 
    % Case 8 ** 
    k.Ta  = 273.15+20;                % [K]
    k.xa      = 2103.3./k.pa; % 2337 * 0.9 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;    
    
elseif tempCond == 25 
    % Case 8 ** 
    k.Ta  = 273.15+25;                % [K]
    k.xa      = 2850.3./k.pa; % 3167 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
    
elseif tempCond == 2520
    % Tamb = 25 C, RH = 20 %
    k.Ta  = 273.15+25;                % [K]
    k.xa      = 633.4./k.pa; % 3167 * 0.2 % 20 % RH condition 
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;

end      


% Deep reindeer's body temperature
if noseGeom == 0
    k.T_body  = 273.15 + 36.3;          %38.4; %***

elseif noseGeom == 2 | noseGeom == 4    % 2, 4: arctic 
    k.T_body  = 273.15 + 36.3; % 35.9; %%  38.4; %ref [***]
    
elseif noseGeom == 1 | noseGeom == 3    % 1, 3: subtropical
    k.T_body  = 273.15 + 35.7; %35.7; %%  38.4; %ref [***]

elseif noseGeom == 5                    %  subtropical from Huntley 1984
    k.T_body  = 273.15 + 35.7; %35.7; %%  38.4; %ref [***]  
    
elseif noseGeom == 6 | noseGeom == 7 | noseGeom == 8  % 6, 7, 8 
    k.T_body  = 273.15 + 36.3; % 35.9;    %%  38.4; %ref [***]
    
elseif noseGeom == 9 | noseGeom == 10 | noseGeom == 11  % 9, 10, 11 
    k.T_body  = 273.15 + 36.3; % 35.9;    %%  38.4; %ref [***]
    
elseif noseGeom == 12  % deep body for subtropical seal 
    k.T_body  = 273.15 + 37; % 35.9;    %%  38.4; %ref [***]

else
    k.T_body  = 273.15 + 37; % for other animals.
    
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
if noseGeom == 0  % reindeer 
    M_seal= 54;           %*** [kg]
    % Respiratory minute volume of air
    V_min     =  0.1; % [l/(min kg_reindeer)]
    k.V         = V_min *1e-3*(M_seal); % [m3 in 1 min]  %***
    M_tot     = k.rhoa*k.V;   % kg in 1 minuto 

elseif noseGeom == 20  % human
    M_seal= 70;           %*** [kg]  !! essential
    % Respiratory minute volume of air
    V_min     = 0.1; % [l/(min kg_reindeer)]
    k.V       = 0.125*60*1e-3; % [m3 in 1 min]  %***
    M_tot     = k.rhoa*k.V;   % kg in 1 minuto %***!! essential

else % seal 
    M_seal= 180;   % reed 1994           %96  r2 ;%54;           %*** [kg]
    % Respiratory minute volume of air
    %V_min     = 0.3967 r2 ; %*** [l/(min kg_seal^0.75)]     0.1; % [l/(min kg_reindeer)]
    %k.V         = V_min *1e-3*(M_seal^0.75)  r2 ; % [m3 in 1 min]  %*** regarding Folkow 1987 and 96 kg, it is 12.16L 
    V_min = 0.3; % folkow 1987
    k.V = rTide*V_min *1e-3*(M_seal)^0.75;
    %k.V = 6.3 * 19.4 % 6.3 L per cycle and 19.4 per min   r2: 12.16L 
    % k.V  = 6.3 r2  % 6.3 L per cycle and 
    M_tot     = k.rhoa*k.V;   % kg in 1 minuto 
end


% Air mass flow
k.F_a     = M_tot/60;       % [kg/s]
% Blood mass flow in veins and arteries (PROPER VALUE HERE: in the whole head is 48e-4!!!!!)
k.F_b   = 2e-04*M_seal/M_reindeer;          % [kg/s] % *** A fraction of the value above, gives velocities that makes sense, see article
% Blood density
k.rho_b = 1000;             % [kg/m3]
k.rho_w = 1000;             % [kg/m3]

% Heat capacities of air, water, and blood (we assume them constant)
k.c_a   = 1003;             % [J/kg K]
k.c_w   = 4186;             % [J/kg K]
k.c_b   = 3620;            % [J/kg K] blood for 3620, 4509 for water 
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
seal = DATA_Csvdata(36, r_shunt);
sealHuntley = DATA_Txtdata(36);


if noseGeom == 0 
    k.L      = k.L_reindeer;
    k.Area_a = k.Area_a_reindeer;
    k.Per_a  = k.Per_a_reindeer;
    
elseif noseGeom == 1 | noseGeom == 12 % 1: subtropical
    k.L      = 1e-3*(max(seal.xSubPosi) - min(seal.xSubPosi));
    k.Area_a = 1e-6*flip(seal.xSubAreaWoshunt);
    k.Per_a  = 1e-3*flip(seal.xSubPeri);

elseif noseGeom == 2 % 2: arctic
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = 1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 3 % 3: subtropical+shunt
    k.L      = 1e-3*(max(seal.xSubPosi) - min(seal.xSubPosi));
    k.Area_a = 1e-6*flip(seal.xSubAreaWshunt);
    k.Per_a  = 1e-3*flip(seal.xSubPeri);
    
elseif noseGeom == 4 % 4: arctic+shunt
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*flip(seal.xArcAreaWshunt); % was Woshunt 
    k.Per_a  = 1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 5 % 5: subtropical from Huntley 1984
    k.L      = 1e-3*(max(sealHuntley.xArcHuntley) - min(sealHuntley.xArcHuntley));
    k.Area_a = 1e-6*flip(sealHuntley.areaArcHuntley);
    k.Per_a  = 1e-3*flip(sealHuntley.periArcHuntley);

elseif noseGeom == 6 % 2: arctic A101 area 1.2
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = (1.2^2)*1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = 1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 7 % 2: arctic A2 peri 1.2
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = (1.2)*1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 8 % 2: arctic A3 area+peri 1.2
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = (1.2^2)*1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = (1.2)*1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 9 % 2: arctic A1 area 0.8 
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = (0.8^2)*1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = 1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 10 % 2: arctic A2 peri 0.8 
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = (0.8)*1e-3*flip(seal.xArcPeri);
    
elseif noseGeom == 11 % 2: arctic A3 area+peri 0.8
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = (0.8^2)*1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = (0.8)*1e-3*flip(seal.xArcPeri); 

elseif noseGeom == 20 % human 
	% all in cm, cm, cm**2
	x_human = [0.02417951923 0.3488976049 0.6888573594 1.013441511 1.341954388 1.684280307 1.995542508 2.335064746 2.654621915 2.978277459 3.311603025 3.652455671 3.968262693 4.295570166 4.626440279 4.949506514 5.281082011 5.935946967 6.269835055 6.585552787 6.937905885 7.258882752 7.573930816 7.909747551 8.255993266].*1e-2;
	p_human = [8.764237274 8.374849188 7.43221642 6.731666311 8.311344957 8.991467436 9.500415957 12.15082473 15.07724606 17.34741344 18.71893933 18.50163987 18.80335031 19.76070491 19.71640874 17.42733818 14.68833709 9.21155446 4.504096415 4.667648122 4.933904204 5.235004856 5.329477201 5.216101676 5.205430526].*1e-2;
	a_human = [1.855691207 1.465998232 1.42139718 1.239247183 1.816360397 1.979607215 2.281317647 2.305186179 2.606896612 2.735603748 3.244247379 3.165106646 3.190194738 3.457060599 3.4818438 3.713865091 2.97908019 3.097725956 2.120553506 2.21472096 2.549446625 2.574839606 2.219991202 2.140850469 2.027474945].*1e-4;

    k.L      = 0.082;               %  8.2 [cm -> m]
    xx_human = 0:k.L/(Ns-1):k.L;
    k.Per_a  = pchip(x_human, p_human, xx_human); %[cm -> m]
    k.Area_a = pchip(x_human, a_human, xx_human); %[cm**2 -> m]
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
%k.d       = 0.7e-03;          % [m]  
k.d       = d_mucus ;             % [m]  % r6
k.Per_w   = k.Per_a;
k.Area_w  = k.d.*k.Per_w;
k.A_w_mean  = trapz(k.x_w,k.Area_w)/k.L;
k.P_w_mean  = trapz(k.x_w,k.Per_w)/k.L;
%k.d_mucus = 5e-4; % r6
k.d_mucus = k.d ;

% Arteries (area and perimeter of the cross sections) [m2 and m]
k.x_art   = linspace(0,k.L,5);
P_art     = pchip(k.x_a,k.Per_a,k.x_art);
A_art     = [5.79 3.88 1.17 1.53 3.85]*1e-5 ;  %***max(k.Area_a)/max(k.Area_a_reindeer)    % area per m of lining 
k.Area_art= A_art.*P_art*factor_rtot;
k.Per_art = [0.772 0.62 0.323 0.33 0.44].*P_art*factor_rtot;
%k.frac_art= [5.3 3.5 1.9 1.5 4.8]*1e-2*factor_rit; % ***
k.frac_art= [5.79 3.88 1.17 1.53 3.85]*1e-2; % ***


k.A_art_mean  = trapz(k.x_art,k.Area_art)/k.L;
k.P_art_mean  = trapz(k.x_art,k.Per_art)/k.L;

% Veins (area and perimeter of the cross sections) [m2 and m]
k.x_ven   = linspace(0,k.L,5);
A_ven     = [6.70 6.08 1.91 2.26 5.50]*1e-4; %***       % area per m of lining 
P_ven     = pchip(k.x_a,k.Per_a,k.x_ven);
k.Area_ven= A_ven.*P_ven*factor_rtot;
k.Per_ven = [3.31 2.76 1.70 1.51 2.20].*P_ven*factor_rtot;
%k.frac_ven= [42 38 25 20 26]*1e-2*factor_rit;***
k.frac_ven= [67.0 60.8 19.1 22.6 55.0]*1e-2;

k.A_ven_mean  = trapz(k.x_ven,k.Area_ven)/k.L;
k.P_ven_mean  = trapz(k.x_ven,k.Per_ven)/k.L;


k.Area_tot= (k.Area_art./k.frac_art);

% Meat 
k.x_m     = linspace(0,k.L,5);
k.dist_art= [0.2 0.11 0.1 0.18 0.64]*1e-03*max(k.Area_a)/max(k.Area_a_reindeer); %*** 
k.dist_ven= [0.4 0.25 0.15 0.15 0.225]*1e-03*max(k.Area_a)/max(k.Area_a_reindeer); %***
k.frac_m  = (1-k.frac_art-k.frac_ven)*factor_rit;
k.Area_m  = k.Area_tot.*k.frac_m;
k.Area_m_mean  = trapz(k.x_m,k.Area_m)/k.L;




%% BREATHING CYCLE
if noseGeom == 0 
    k.fa      = 10;
elseif noseGeom == 20
	k.fa      = 14.3; % Lindemann 14 breaths per min. 
else  
    k.fa      = 19.4*rFreq; %19.4; % Reed and butler 1993   %10.6;   % breaths per min from [***] 
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
% Minimum velocity for 
k.velLimit = vLimit; 

%% factor for increase or decrease resistance at interface. 
k.factor_Ram = factorRam ; % 1 ; % 10;
k.factor_Rij = factorRij ; % 1 ; % 0.1;
k.factor_fric = factorFric; 


k.h_a = heatTranCoef; 
%  d_mucus, factorRam, factorRij


