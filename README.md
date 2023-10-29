# nose-calculations
This repository is about a numerical calculation code package to study hydrodynamic and thermodynamic properties during respiration in various animal noses. 

## Getting started
#### The easiest and simplest way to download all codes for modelling an animal nose is by 'nose-calculations' and run in Matlab. 


1. git clone https://github.com/hyejeonc/nose-calculations.git
2. Open 'Nose_t.m' file in Matlab.
3. Modify 'Nose_t.m' file, i.e., changing the ambient air temperature, you should modify \\\k = Parameters_animals(.. , .. , .. , .. , ... )
```matlab
function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
% Important parameters are tempCond, noseGeom.
% tempCond = the ambient air temperature
% noseGeom = selection of nose geometry regarding each animal. Selected nose geometry should be added in Parameters_animals.
% Explanations of other parameters are described later.


```
4. Depending on the parameters for each animal, modify 'Parameters_animals.m' file.#### Nose_t.m. 
```matlab
% 
elseif tempCond == -30
    % Case 1 
    k.Ta  = 273.15-30;                % [K]
    k.xa      = 34.2/k.pa; % 38 * 0.9
    Ma        = k.xa.*ktp1.Mw(1)+0.79*(1-k.xa).*ktp1.Mw(2)+0.21*(1-k.xa).*ktp1.Mw(3);
    k.wa      = k.xa*ktp1.Mw(1)/Ma;

% add 
elseif noseGeom == 2 % 2: arctic
    k.L      = 1e-3*(max(seal.xArcPosi) - min(seal.xArcPosi));
    k.Area_a = 1e-6*flip(seal.xArcAreaWoshunt);
    k.Per_a  = 1e-3*flip(seal.xArcPeri);
    
```

5. 



## A main file for running a simulation. 

#### Parameters.m
## To change 'default' physical and thermodynamic parameters, you should change it, i.e., heat capacity of mucus lining. 
