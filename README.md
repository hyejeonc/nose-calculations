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
The output file '***.mat' consists of a struct named Nose, and this has the following information. 
```matlab
        Nose.Ta = array, temperature of air subsystem.
        Nose.Tw = np.array(Nose['Tw'])
        Nose.Tit = np.array(Nose['Tm'])
        Nose.Tm = np.array(Nose['Tm'])
        Nose.Tven = np.array(Nose['Tven'])
        Nose.wa = np.array(Nose['wa'])
        Nose.x_a = np.array(Nose['x_a'])
        Nose.Tart_art = np.array(Nose['Tart_art'])
        Nose.Tart = np.array(Nose['Tart'])
        Nose.M_tot = np.array(Nose['M_tot'])
        Nose.rho = np.array(Nose['rho'])
        Nose.rho_vap = np.array(Nose['rho_vap'])
        Nose.rho_dry = np.array(Nose['rho_dry'])
        Nose.va = np.array(Nose['va'])
        Nose.F_vap = np.array(Nose['F_vap'])
        Nose.Fa = np.array(Nose['Fa'])
        Nose.vart = np.array(Nose['vart'])
        Nose.vven = np.array(Nose['vven'])
        Nose.R_aw = np.array(Nose['R_aw'])
        Nose.R_wm = np.array(Nose['R_wm'])
        Nose.R_art = np.array(Nose['R_art'])
        Nose.R_ven = np.array(Nose['R_ven'])
        Nose.R_mu = np.array(Nose['R_mu'])
        Nose.R_qmu = np.array(Nose['R_qmu'])
        Nose.h_a = np.array(Nose['h_a'])
        Nose.h_wa = np.array(Nose['h_wa'])
        Nose.h_w = np.array(Nose['h_w'])
        Nose.X_T = np.array(Nose['X_T'])
        Nose.X_mu = np.array(Nose['X_mu'])
        Nose.Jw = np.array(Nose['Jw'])
        Nose.Jq_a = np.array(Nose['Jq_a'])
        Nose.Jq_w = np.array(Nose['Jq_w'])
        Nose.Jq_art = np.array(Nose['Jq_art'])
        Nose.Jq_ven = np.array(Nose['Jq_ven'])
```


## Trouble shooting



#### Parameters.m
## To change 'default' physical and thermodynamic parameters, you should change it, i.e., heat capacity of mucus lining. 
