# nose-calculations
This repository is about a numerical calculation code package to study hydrodynamic and thermodynamic properties during respiration in various animal noses. 

## Getting started
#### The easiest and simplest way to download all codes for modelling an animal nose is by 'nose-calculations' and run in Matlab. 


1. git clone https://github.com/hyejeonc/nose-calculations.git
2. open 'Nose_t.m' file in Matlab.
3. modify 'Nose_t.m' file, i.e., changing the ambient air temperature, you should modify `k = Parameters_animals(.. , .. , .. , .. , ... )`
   
```matlab
k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
% Important parameters are tempCond, noseGeom.
% tempCond = the ambient air temperature [degree Celsius], for example, -30.
% noseGeom = selection of nose geometry regarding each animal. Selected nose geometry should be added in 'Parameters_animals.m'.
```

4. Depending on the parameters for each animal, modify 'Parameters_animals.m' file.
   To add an additional ambient temperature condition, add elseif paragraph here:
   k.xa = water mass fraction 
    which refers from the saturate water vapour in the air from https://www.engineeringtoolbox.com/moist-air-properties-d_1256.html
    You can modify relative humidity by changing the multiplying factor, 0.9 (RH = 90 %). 

```matlab
elseif tempCond == -30
    % Case 1 
    k.Ta  = 273.15-30;                % [K]
    k.xa      = 34.2/k.pa; % 38 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;
```

To add geometrical information, you should add a length, cross-sectional area and perimeter information here: 
add another elseif paragraph to avoid a confusion. 

```matlab
elseif noseGeom == 2 % 2: arctic
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi)); [m]
    k.Area_a = 1e-6*flip(seal.xArcAreaWoshunt); [m^2]
    k.Per_a  = 1e-3*flip(seal.xArcPeri); [m]
```

5. Check the filename to save 'saveName.mat' file in 'Nose_t.m' file.

```matlab
k.fileName = ['saveName','.mat'];
```

6. Run 'Nose_t.m' in Matlab.
7. The output file 'saveName.mat' file is saved automatically after one iteration in the same directory of 'Nose_t.m'.



## Details of the main code - 'Nose_t.m'
```matlab
function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
% Important parameters are tempCond, noseGeom.
% tempCond = the ambient air temperature [degree Celsius], for example, -30.
% noseGeom = selection of nose geometry regarding each animal. Selected nose geometry should be added in 'Parameters_animals.m'.
```


## Data saved in the output file 
The output file '***.mat' consists of a struct named 'Nose' and this has the following information. 
All is saved as [number of iteration * 2 * Nt, Ns] array. 
For example, the solution (in the last iteration) of the temperature of the air subsystem during inhalation is: Nose.Ta(end-2*Nt:end-Nt,:)
For exhalation: Nose.Ta(end-Nt:,:)

```matlab
        Nose.Ta = temperature of air subsystem.
        Nose.Tw = temperature of mucus subsystem.
        Nose.Tit = temperature of interstitial tissue subsystem.
        Nose.Tven = temperature of vein subsystem.
        Nose.wa = water mass fraction in air subsystem.
        Nose.x_a =  water mol fraction in air subsystem.
        Nose.Tart = temperature of artery subsystem.
        Nose.M_tot = total mass of air 
        Nose.rho = density of air 
        Nose.rho_vap = density of water vapour in the air
        Nose.rho_dry = density of dry air (density after excluding water vapour)
        Nose.va = flow velocity of air
        Nose.F_vap = mass flow rate of dry air
        Nose.Fa = mass flow rate of air 
        Nose.R_aw = interface resistivity air-mucus (heat, Rqq)
        Nose.R_wm = interface resistivity mucus-interstitial tissue (heat)
        Nose.R_art = interface resistivity interstitial tissue-artery (heat)
        Nose.R_ven = interface resistivity interstitial tissue-vein (heat)
        Nose.R_mu = interface resistivity of water mass air-mucus (mass)
        Nose.R_qmu = coupling coefficient (heat-mass)
        Nose.h_a = enthalpy of air 
        Nose.h_wa = enthalpy of water vapour
        Nose.h_w = enthalpy of liquid water
        Nose.X_T = driving force of heat transport (difference between reciprocal of temperature)
        Nose.X_mu = driving force of mass transport (chemical potential difference)
        Nose.Jw = water mass flux air-mucus subsystem 
        Nose.Jq_a = heat flux of air-mucus subsystem
        Nose.Jq_w = heat flux of mucus-interstitial tissue subsystem
        Nose.Jq_art = heat flux of interstitial tissue - artery subsystem
        Nose.Jq_ven = heat flux of interstitial tissue - vein subsystem
```

## Trouble shooting



#### Parameters.m
## To change 'default' physical and thermodynamic parameters, you should change it, i.e., heat capacity of mucus lining. 
