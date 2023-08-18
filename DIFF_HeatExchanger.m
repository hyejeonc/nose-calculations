function der=DIFF_HeatExchanger(t,y)

global k ktp1 ktp2
y=y';
t    = t.*k.Scales.t;
Ta   = y(1:k.Ns).*k.Scales.T;
Tm   = y(k.Ns+1:2*k.Ns).*k.Scales.T;
Tit  = y(2*k.Ns+1:3*k.Ns).*k.Scales.T;
Tart = y(5*k.Ns+1:5*k.Ns+k.Ns_art).*k.Scales.T;
Tven = y(3*k.Ns+1:4*k.Ns).*k.Scales.T;

% Water mass fraction in air
w_vap= y(4*k.Ns+1:5*k.Ns);      
M_dry= 0.79*ktp1.Mw(2) + 0.21*ktp1.Mw(3);
M_vap= ktp1.Mw(1);
x_vap= w_vap.*M_dry./(M_vap.*(1-w_vap)+w_vap.*M_dry);
Ma   = x_vap.*M_vap + (1-x_vap).*M_dry;
xa   = [x_vap; (1-x_vap)*0.79; (1-x_vap)*0.21];

for ii = 1:k.Ns
    [Zfac_struct,~] = TP_Zfac_a(Ta(ii),k.pa,xa(:,ii),2);
    rho_a(1,ii)     = Zfac_struct.c.*Ma(ii);
    rho_vap(1,ii)   = rho_a(ii).*w_vap(ii);
    rho_dry(1,ii)   = rho_a(ii)-rho_vap(1,ii);
end

% Heat capacity of air
c_a  = zeros(size(Ta));
for ii = 1:k.Ns
    cp0  = IG_Cp_a(Ta(ii),xa(:,ii));
    cp_res      = CB_Enthalpy_a(Ta(ii),k.pa,xa(:,ii));
    c_a(1,ii)   = (cp0 + cp_res.dT)./Ma(ii);
end

% Heat capacity of mucus
c_w  = zeros(size(Tm));
for ii = 1:k.Ns
    cp0         = IG_Cp_w(Tm(ii));
    cp_res      = CB_Enthalpy_w(Tm(ii),k.pa);
    c_w(1,ii)   = (cp0 + cp_res.dT)./ktp2.Mw;
end

% Chemical potential & enthalpies 
for ii  = 1:k.Ns
    mu_str      = TP_ChemicalPotential_a(Ta(ii),k.pa,xa(:,ii),1);
    mu_a(1,ii)  = mu_str.mu(1);
    h_id        = IG_H_a(Ta(ii),xa(:,ii));
    h_res       = CB_Enthalpy_a(Ta(ii),k.pa,xa(:,ii));
    h_tot(1,ii) = (h_id.H + h_res.H)./Ma(ii);
    h_i(:,ii)   = (h_id.dZ + h_res.dN)./ktp1.Mw;
    h_a(1,ii)   = h_i(1,ii);
 
    mu_str      = TP_ChemicalPotential_w(Tm(ii),k.pa,1);
    mu_m(1,ii)  = mu_str.mu;
    h_id        = IG_H_w(Tm(ii));
    h_res       = CB_Enthalpy_w(Tm(ii),k.pa);
    h_w(1,ii)   = (h_id.H + h_res.H)./ktp2.Mw;
end
 
% Geometry parameter interpolation
gamma_a  = k.P_a_mean*ones(size(k.zz));
gamma_art= k.P_art_mean*ones(size(k.zz_art));
gamma_art_int= k.P_art_mean*ones(size(k.zz));
gamma_ven= k.P_ven_mean*ones(size(k.zz));
gamma_ven_a= k.P_ven_mean*ones(size(k.zz_art));
A_a      = k.A_a_mean*ones(size(k.zz));    % [m2]
A_it     = k.Area_m_mean*ones(size(k.zz));
A_m      = k.A_w_mean*ones(size(k.zz));
A_art    = k.A_art_mean*ones(size(k.zz_art));     % [m2]
A_art_int= k.A_art_mean*ones(size(k.zz));     % [m2]
A_ven    = k.A_ven_mean*ones(size(k.zz));     % [m2]
A_ven_a  = k.A_ven_mean*ones(size(k.zz_art));     % [m2]

% Hydraulic diameters
D_a      = 4.*A_a./gamma_a;
D_art    = 4.*A_art_int./gamma_art_int;
D_ven    = 4.*A_ven./gamma_ven;
D_art_a  = 4.*A_art./gamma_art;
D_ven_a  = 4.*A_ven_a./gamma_ven_a;

% Air flow
F_dry    = k.A_breathing*sin(pi/k.T_breathing*t).*ones(1,k.Ns);
v        = F_dry./(rho_dry.*A_a);
F_vap    = A_a.*v.*rho_vap;

% Temperature derivatives
dTadz   = dss020(k.zz(1),k.zz(k.Ns),k.Ns,Ta  , 1);
dTartdz = dss020(k.zz(1),k.zz(k.Ns),k.Ns_art,Tart,-1);
dTvendz = dss020(k.zz(1),k.zz(k.Ns),k.Ns,Tven,1);
dFvapdz = dss020(k.zz(1),k.zz(k.Ns),k.Ns,F_vap,1);
dwvapdz = dss020(k.zz(1),k.zz(k.Ns),k.Ns,w_vap,1);
dw2dz   = -0.79.*ktp1.Mw(2)./M_dry.*dwvapdz;
dw3dz   = -0.21.*ktp1.Mw(3)./M_dry.*dwvapdz;
dwadz   = [dwvapdz; dw2dz; dw3dz];

Tart_int = pchip(k.zz_art, Tart, k.zz);
Tit_art = pchip(k.zz, Tit, k.zz_art);
% Fluxes between subsystems
for ii = 1:k.Ns
    R_am        = CALC_R_a_m(Ta(ii),xa(:,ii),Tm(ii),D_a(ii),v(ii));
    X           = [1./Tm(ii)-1./Ta(ii);
                 -(mu_m(ii)./Tm(ii)-mu_a(ii)./Ta(ii))+h_a(1,ii).*(1./Tm(ii)-1./Ta(ii))];
    JJ         = R_am\X;
    Jw(1,ii)   = JJ(2);
    Jq_a(1,ii) = JJ(1);
    % Fluxes between subsystems
    R_ij           = CALC_R_ij(D_art(ii),D_ven(ii));
    Jq_w(1,ii)     = 1./R_ij.m_it   .* (1./Tit(ii)-1./Tm(ii));
    Jq_art(1,ii)   = 1./R_ij.it_art .* (1./Tart_int(ii)-1./Tit(ii));
    Jq_ven(1,ii)   = 1./R_ij.it_ven .* (1./Tven(ii)-1./Tit(ii));
end

for ii = 1:k.Ns_art
    % Fluxes between subsystems
    R_ij           = CALC_R_ij(D_art_a(ii),D_ven_a(ii));
    Jq_art_a = 1./R_ij.it_art .* (1./Tart-1./Tit_art);
end
    
drhovapdt = (-dFvapdz./A_a - gamma_a./A_a.*Jw);
dwvapdt   = drhovapdt.*(1-w_vap)./rho_a;
dw2dt     = -0.79.*ktp1.Mw(2)./M_dry.*dwvapdt;
dw3dt     = -0.21.*ktp1.Mw(3)./M_dry.*dwvapdt;
dwadt     = [dwvapdt; dw2dt; dw3dt];

for ii = 1:k.Ns
    term_w(1,ii) = sum(h_i(:,ii).*dwadt(:,ii),1) + v(:,ii).*sum(h_i(:,ii).*dwadz(:,ii),1);
end

dTadt   = (-v.*dTadz - gamma_a./(A_a.*rho_a.*c_a).*(Jq_a-Jw.*(h_tot-h_i(1,:))) - 1./c_a.*(term_w) )./k.Scales.T*k.Scales.t;
% dTmdt   = ( gamma_a .*(Jq_a- Jq_w+ Jw.*(0))./(k.rho_w.*c_w.*A_m)                               )./k.Scales.T*k.Scales.t;
dTmdt   = ( gamma_a .*(Jq_a- Jq_w+ Jw.*(h_i(1,:)-h_w(1,:)))./(k.rho_w.*c_w.*A_m)                               )./k.Scales.T*k.Scales.t;
dTitdt  = ((gamma_a .*Jq_w -gamma_art_int.*Jq_art-gamma_ven.*Jq_ven)./(k.rho_b.*k.c_b.*A_it))./k.Scales.T*k.Scales.t;
dTartdt = ( k.F_b./(A_art.*k.rho_b).*dTartdz + gamma_art.*Jq_art_a./(A_art.*k.rho_b.*k.c_b))./k.Scales.T*k.Scales.t;
dTvendt = (-k.F_b./(A_ven.*k.rho_b).*dTvendz + gamma_ven.*Jq_ven./(A_ven.*k.rho_b.*k.c_b))./k.Scales.T*k.Scales.t;

dwvapdt   = dwvapdt.*k.Scales.t;
w_ref     = k.wa;

% Boundary conditions in z = 0
res_wa    = w_vap(1) - w_ref;
res_Ta1   = (Ta(1)-k.Ta)./k.Scales.T;
res_Tven1 = (Tven(1)-Tart(1))./k.Scales.T;

% Boundary conditions in z = L
res_Tm2   = (Tm(end)-k.T_body)  ./k.Scales.T;
res_Tit2  = (Tit(end)-k.T_body) ./k.Scales.T;
res_Tart2 = (Tart(end)-k.T_body)./k.Scales.T;

der       = [  res_Ta1 dTadt(2:end) dTmdt(1:end-1) res_Tm2 dTitdt(1:end-1) res_Tit2  res_Tven1 dTvendt(2:end) res_wa dwvapdt(2:end) dTartdt(1:end-1) res_Tart2]';