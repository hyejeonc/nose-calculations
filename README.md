# nose-calculations
This repository is about a numerical calculation code package to study hydrodynamic and thermodynamic properties during respiration in various animal noses. 

## Getting started
#### The easiest and simplest way to download all codes for modelling an animal nose is by 'nose-calculations' and run in Matlab. 


1. git clone https://github.com/hyejeonc/nose-calculations.git
2. Open 'Nose_t.m' file in Matlab.
3. Modify 'Nose_t.m' file, i.e., changing the ambient air temperature, you should modify `k = Parameters_animals(.. , .. , .. , .. , ... )`
   
```matlab
function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
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
please add another elseif paragraph to avoid a confusion. 

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

6. Run 'Nose_t.m'.
7. The output file 'saveName.mat' file is save in the same directory of 'Nose_t.m'.



## Details of the main code - 'Nose_t.m'
```matlab
function k = Parameters_animals(tempCond, noseGeom, Ns, Ns_art, N_cycle, r_shunt, d_mucus, factorRam, factorRij, factorFric, heatTranCoef, factor_rit, factor_rtot, rTide, rFreq, velLim)
% Important parameters are tempCond, noseGeom.
% tempCond = the ambient air temperature [degree Celsius], for example, -30.
% noseGeom = selection of nose geometry regarding each animal. Selected nose geometry should be added in 'Parameters_animals.m'.
```


## Data saved in the output file 
The output file '***.mat' is consisted of a struct named Nose, and this have the following information. 
```matlab
       self.Ta = np.array(Nose['Ta'])
        self.Tw = np.array(Nose['Tw'])
        self.Tit = np.array(Nose['Tm'])
        self.Tm = np.array(Nose['Tm'])
        self.Tven = np.array(Nose['Tven'])
        self.wa = np.array(Nose['wa'])
        self.x_a = np.array(Nose['x_a'])
        self.Tart_art = np.array(Nose['Tart_art'])
        self.Tart = np.array(Nose['Tart'])
        self.M_tot = np.array(Nose['M_tot'])

        self.rho = np.array(Nose['rho'])
        self.rho_vap = np.array(Nose['rho_vap'])

        self.rho_dry = np.array(Nose['rho_dry'])
        self.va = np.array(Nose['va'])
        self.F_vap = np.array(Nose['F_vap'])
        self.Fa = np.array(Nose['Fa'])
        self.vart = np.array(Nose['vart'])
        self.vven = np.array(Nose['vven'])
        self.R_aw = np.array(Nose['R_aw'])
        self.R_wm = np.array(Nose['R_wm'])
        self.R_art = np.array(Nose['R_art'])
        self.R_ven = np.array(Nose['R_ven'])
        self.R_mu = np.array(Nose['R_mu'])
        self.R_qmu = np.array(Nose['R_qmu'])
        self.h_a = np.array(Nose['h_a'])
        self.h_wa = np.array(Nose['h_wa'])
        self.h_w = np.array(Nose['h_w'])
        self.X_T = np.array(Nose['X_T'])
        self.X_mu = np.array(Nose['X_mu'])
        self.Jw = np.array(Nose['Jw'])
        self.Jq_a = np.array(Nose['Jq_a'])
        self.Jq_w = np.array(Nose['Jq_w'])
        self.Jq_art = np.array(Nose['Jq_art'])
        self.Jq_ven = np.array(Nose['Jq_ven'])
        self.sigma = np.array(Nose['sigma'])
```


## Trouble shooting



#### Parameters.m
## To change 'default' physical and thermodynamic parameters, you should change it, i.e., heat capacity of mucus lining. 
