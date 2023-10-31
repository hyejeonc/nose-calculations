function k = ScalesParameters(k)
 

% Scaling factors of the variables:
k.Scales.T = k.T_body;          % [K]
k.Scales.p = k.pa;              % [Pa]
k.Scales.rho = k.rhoa;          % [m]
k.Scales.F = k.F_a.*1e-2;       % [m]
k.Scales.L = k.L;               % [kmol/m3]
k.Scales.t = k.T_breathing;     % [kmol/m3]