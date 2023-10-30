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
k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
% Important parameters are tempCond, noseGeom.
% tempCond = the ambient air temperature [degree Celsius], for example, -30.
% noseGeom = selection of nose geometry regarding each animal. Selected nose geometry should be added in 'Parameters_animals.m'.
```

## Data saved in the output file 
The output file 'saveName.mat' consists of a struct named 'Nose' and this has the following information. 

All is saved as [number of iteration * 2 * Nt, Ns] array. 

For example, the solution (in the last iteration) of the temperature of the air subsystem during inhalation is: Nose.Ta(end-2*Nt:end-Nt,:)

For exhalation: Nose.Ta(end-Nt:,:)

```matlab
Nose.Ta      = temperature of air subsystem.
Nose.Tw      = temperature of mucus subsystem.
Nose.Tit     = temperature of interstitial tissue subsystem.
Nose.Tven    = temperature of vein subsystem.
Nose.wa      = water mass fraction in air subsystem.
Nose.x_a     = water mol fraction in air subsystem.
Nose.Tart    = temperature of artery subsystem.
Nose.M_tot   = total mass of air 
Nose.rho     = density of air 
Nose.rho_vap = density of water vapour in the air
Nose.rho_dry = density of dry air (density after excluding water vapour)
Nose.va      = flow velocity of air
Nose.F_vap   = mass flow rate of dry air
Nose.Fa      = mass flow rate of air 
Nose.R_aw    = interface resistivity air-mucus (heat, Rqq)
Nose.R_wm    = interface resistivity mucus-interstitial tissue (heat)
Nose.R_art   = interface resistivity interstitial tissue-artery (heat)
Nose.R_ven   = interface resistivity interstitial tissue-vein (heat)
Nose.R_mu    = interface resistivity of water mass air-mucus (mass)
Nose.R_qmu   = coupling coefficient (heat-mass)
Nose.h_a     = enthalpy of air 
Nose.h_wa    = enthalpy of water vapour
Nose.h_w     = enthalpy of liquid water
Nose.X_T     = driving force of heat transport (difference between reciprocal of temperature)
Nose.X_mu    = driving force of mass transport (chemical potential difference)
Nose.Jw      = water mass flux air-mucus subsystem 
Nose.Jq_a    = heat flux of air-mucus subsystem
Nose.Jq_w    = heat flux of mucus-interstitial tissue subsystem
Nose.Jq_art  = heat flux of interstitial tissue - artery subsystem
Nose.Jq_ven  = heat flux of interstitial tissue - vein subsystem
```


## Change simulation settings
Most of simulation settings can be changed from `k = Parameters_animals(.. , .. , .. , .. , ... )` and 'Parameters_animals.m'.

```matlab
k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
% Ns           = discretization number in space 
% Ns_art       = discretization number in space for artery only
% N_cycle      = maximum number of cycle
% r_shunt      = a fraction to segmentize nose length. default = 0.5. This means you segmentize the nose length where the volume ratio between maxilloturbinate and total air pathway in a nose (maxilloturbinate+olfactory path) is 0.5.
% d_mucus      = thickness of mucus lining [m]
% factorRam    = multiplying factor for resistivity between air-mucus
% factorRij    = multiplying factor for resistivity between two subsystems which does not relate to mass transport (i.e., interstitial tissue-artery, interstitial tissue-vein, mucus-interstitial tissue)
% factorFric   = multiplying factor for the friction of the air in a turbinate 
% heatTranCoef = constant heat transfer coefficient regardless of nose geometry 
% factor_rit   = multiplying factor for volume ratio of interstitial tissue; factor_rit * interstitial tissue/(interstitial tissue + artery + vein)
% factor_rtot  = multiplying factor total volume of (interstitial tissue + artery + vein)
% rTide        = multiplying factor for tidal volume
% rFreq        = multiplying factor for respiration frequency
% velLim       = the minima of air velocity, which avoids Reynolds number diverging to infinity
```


## Troubleshooting
1. The solution does not converge.
   - [ ] Increase Ns and/or Ns_art in k = Parameters_animals()
   - [ ] Check residuals of several variables in Nose struct, and increase N_cycle in `k = Parameters_animals()`
2. Message from Matlab 'To run this file, you can either change the MATLAB current folder or add its folder to the MATLAB path'.
   - [ ] Be sure that the directory for 'Nose_t.m' is with 'AIR_IG' and 'Water_EoS' folders. 


## Citation 
Add the following citations when you use this repo for scientific purposes. 

[1] Magnanelli, Elisa, et al. "The nasal geometry of the reindeer gives energy-efficient respiration." Journal of Non-Equilibrium Thermodynamics 42.1 (2017): 59-78.

[2] Solberg, Simon Birger Byremo, et al. "Energy efficiency of respiration in mature and newborn reindeer." Journal of Comparative Physiology B 190.4 (2020): 509-520.

[3] Cheon, Hyejeong Lee, et al. "Structure-function relationships in the nasal cavity of Arctic and subtropical seals." Biophysical Journal (2023): XXX-XXX.
