function h_a = CALC_R_a_m(T_a, xa,T_m,D_a,v_a)


%%
%% r1 = k.factor_Ram * Rqq
%% r2 = k.factor_Ram * all R elements 
%% r3 = change from Turbulent to Laminar 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_a = abs(v_a);
global k ktp1
RR = load('Coefficients.txt');

if v_a<=1e-4
    v_a=1e-4;
end

Ma   = sum(xa.*ktp1.Mw);
R           = ktp1.R./ktp1.Mw(1);             % J/kgK
[Zfac,~] = TP_Zfac_a(T_a,k.pa,xa,2);
rho_a    = Zfac.c*Ma;

% Dynamic viscosities (Sutherland's law for air)    [kg/(m s)]   
mu         = 1.458e-6*T_a^(3/2)/(T_a+110.4); % for air 
mu_wall    = 1.458e-6*T_m^(3/2)/(T_m+110.4); % for mucus wall 
% Kynematic Viscosity [m2/s]
nu_a= mu./rho_a;
Re_a= v_a *D_a/nu_a;



% disp("====================")
% disp("====================")
% disp("====================")
% % disp("====================")
% disp("Re_a number is == ")
% disp(Re_a)
% disp("====================")
% disp("v_a number is == ")
% disp(v_a)
% disp("====================")
% disp("D_a number is == ")
% disp(D_a)
% disp("====================")
% disp("nu_a number is == ")
% disp(nu_a)
% disp("====================")
% disp("====================")
% disp("====================")
% disp("====================")

Pr_a= nu_a*k.rhoa*k.c_a/k.k_a;
f = 47.78/Re_a*(1+0.127*Re_a^(0.489)) * k.factor_fric;
% Nu_a= 0.023*Re_a^(0.8)*Pr_a^(1/3)
% Nu_a= 0.125*f*Re_a*Pr_a^(1/3);
L = k.L;
Nu_a = 1.86 *(Re_a*Pr_a*D_a/L)^(1/3)*(mu/mu_wall);


h_a = Nu_a*k.k_a./D_a;  % J/m^2sK
%h_a = 200;



v        = 12.7;
M        = ktp1.Mw(1);
D_AB     = T_a^1.75*1e-3/(k.pa/101325*(v^(1/3)+2.73)^2)*(1/M+0.0345)^(1/2)*1e-4;
Sc       = nu_a/D_AB;
% Sider-Tate correlation for mass convection into pipes (it should be 0.023)
Sh       = 0.023*Re_a^(0.8)*Sc^(1/3)/5;
% Sh       = 0.125*f*Re_a*Sc^(1/3);
k_c      = D_AB./D_a.*Sh;
r_22 = R/(k_c.*rho_a);

TT  = RR(:,1);

disp('LAMINAR === ')

disp('1/(h_a*T_a^2) == ')
disp(1/(h_a*T_a^2))

disp("interp1(TT,RR(:,2),T_m,'pchip','extrap') == ")
disp(interp1(TT,RR(:,2),T_m,'pchip','extrap'))

R11 = interp1(TT,RR(:,2),T_m,'pchip','extrap')+ (1/(h_a*T_a^2));
R12 = interp1(TT,RR(:,3),T_m,'pchip','extrap');
R22 = interp1(TT,RR(:,4),T_m,'pchip','extrap') + r_22;

% disp('Ta == ')
% disp(T_a)
% 
% disp('(1/(h_a*T_a^2)')
% disp(1/(h_a*T_a^2))
% 
% disp('h_a')
% disp(h_a)
% 
% disp('R11')
% disp(R11)
% 
disp("R12 = interp1(TT,RR(:,3),T_m,'pchip','extrap');")
disp(interp1(TT,RR(:,3),T_m,'pchip','extrap'))

disp("R22 = interp1(TT,RR(:,4),T_m,'pchip','extrap') + r_22")
disp(interp1(TT,RR(:,4),T_m,'pchip','extrap') + r_22)


%% r1 
% R   = [ R11 R12 ;
R   = [ k.factor_Ram*R11 k.factor_Ram*R12 ; % r2
        k.factor_Ram*R12 k.factor_Ram*R22]; % r2

    

