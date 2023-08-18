function Sol = CALC_OnlyAir
% r0 = r1 : correct dAvendz, dTvendz

global k ktp1 ktp2

k = ScalesParameters(k);

wa      = k.wa;
wb      = k.wb;

if k.initCond == 0
    % Initial conditions
    disp("Initialize temperature and iteration starts! ")
    T0a     = linspace(k.Ta,k.T_body,k.Ns)./k.Scales.T;
    T0w     = linspace(k.T_body,k.T_body,k.Ns)./k.Scales.T;
    T0m     = linspace(k.T_body,k.T_body,k.Ns)./k.Scales.T;
    T0art   = linspace(k.T_body,k.T_body,k.Ns_art)./k.Scales.T;
    T0ven   = linspace(k.T_body,k.T_body,k.Ns)./k.Scales.T;
    w0a     = linspace(wa,wb,k.Ns);
    y0  = [ T0a T0w T0m  T0ven w0a T0art]; 
    
else 
    load(k.initFileName);
    disp("Load file and iteration starts! ")
    T0a     = smooth(Nose.Ta(end,:))./k.Scales.T;
    T0w     = smooth(Nose.Tw(end,:))./k.Scales.T;
    T0m     = smooth(Nose.Tm(end,:))./k.Scales.T;
    T0art   = pchip(k.zz,smooth(Nose.Tart(end,:)),k.zz_art)./k.Scales.T;
    T0ven   = smooth(Nose.Tven(end,:))./k.Scales.T;
    w0a     = smooth(Nose.wa(end,:));
    y0  = [ T0a' T0w' T0m'  T0ven' w0a' T0art];


end


%load('Ref_0.mat');
%T0a     = smooth(Ref_Case.Ta(end,:))./k.Scales.T;
%T0w     = smooth(Ref_Case.Tw(end,:))./k.Scales.T;
%T0m     = smooth(Ref_Case.Tm(end,:))./k.Scales.T;
%T0art   = pchip(k.zz,smooth(Ref_Case.Tart(end,:)),k.zz_art)./k.Scales.T;
%T0ven   = smooth(Ref_Case.Tven(end,:))./k.Scales.T;
%w0a     = smooth(Ref_Case.wa(end,:));
%y0  = [ T0a' T0w' T0m'  T0ven' w0a' T0art];



for i = 1 : k.N_cycle
    n = i-1
    %     Time discretization
    M     = eye(5*k.Ns+ k.Ns_art,5*k.Ns+ k.Ns_art);
    M(1,1)= 0;
    M(k.Ns*3+1,k.Ns*3+1)= 0;
    M(k.Ns*4+1,k.Ns*4+1)= 0;
    M(k.Ns*2,k.Ns*2)= 0;
    M(k.Ns*3,k.Ns*3)= 0;
    M(k.Ns*5+k.Ns_art,k.Ns*5+k.Ns_art)= 0;
    Options = odeset('Mass',M);
    
    tt      = linspace(2.*n.*k.T_breathing,2.*n.*k.T_breathing+k.T_breathing,k.Nt);
    
    [t, x] = ode15s(@DIFF_OnlyAir,tt./k.Scales.t,y0,Options);
    
    disp("IN == ")
    disp(t)
    disp(x)
    % Time
    Sol.t(2*n*k.Nt+1:2*n*k.Nt+k.Nt)   = t*k.Scales.t;
    Sol.x   = k.zz;
    
    % Tempeartures
    Sol.Ta(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)  = x(:,1:k.Ns).*k.Scales.T;
    Sol.Tw(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)  =(x(:,k.Ns+1:2*k.Ns)).*k.Scales.T;
    Sol.Tm(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)  = x(:,2*k.Ns+1:3*k.Ns).*k.Scales.T;
    Sol.Tven(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)= x(:,3*k.Ns+1:4*k.Ns).*k.Scales.T;
    Sol.wa(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)  = x(:,4*k.Ns+1:5*k.Ns);
    Sol.x_a(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)  = Sol.wa(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:).*k.M_dry./(ktp1.Mw(1).*(1-Sol.wa(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:))+Sol.wa(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:).*k.M_dry);
    Sol.Tart_art(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)= x(:,5*k.Ns+1:5*k.Ns+k.Ns_art).*k.Scales.T;
    Sol.Tart(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:)= pchip(k.zz_art,Sol.Tart_art(2*n*k.Nt+1:2*n*k.Nt+k.Nt,:),Sol.x);
    
    gamma_a  = pchip(k.x_a,k.Per_a,k.zz);
    gamma_art= pchip(k.x_art,k.Per_art,k.zz);
    gamma_ven= pchip(k.x_ven,k.Per_ven,k.zz);

    A_a      = pchip(k.x_a,k.Area_a,k.zz);
    A_w      = pchip(k.x_w,k.Area_w,k.zz);
    A_m      = pchip(k.x_m,k.Area_m,k.zz);
    A_art    = pchip(k.x_art,k.Area_art,k.zz);
    A_ven    = pchip(k.x_ven,k.Area_ven,k.zz);
    
    dAadz    = dss020(k.zz(1),k.zz(k.Ns),k.Ns,A_a  , 1);
    dAwdz    = dss020(k.zz(1),k.zz(k.Ns),k.Ns,A_w  , 1);
    dAmdz    = dss020(k.zz(1),k.zz(k.Ns),k.Ns,A_m  , 1);
    dAartdz  = dss020(k.zz(1),k.zz(k.Ns),k.Ns,A_art  , 1);
    dAvendz  = dss020(k.zz(1),k.zz(k.Ns),k.Ns,A_ven  , 1);
    
    cos_a    = 1.*cos(atan(dAadz));
    cos_w    = 1.*cos(atan(dAwdz));
    cos_m    = 1.*cos(atan(dAmdz));
    cos_art  = 1.*cos(atan(dAartdz));
    cos_ven  = 1.*cos(atan(dAvendz));
    
    for jj = 1:k.Nt
        for ii= 1:k.Ns
            xa          = [Sol.x_a(jj,ii); (1-Sol.x_a(jj,ii))*0.79; (1-Sol.x_a(jj,ii))*0.21];
            [Zfac_struct,~] = TP_Zfac_a(Sol.Ta(jj,ii),k.pa,xa,2);
            Sol.M_tot(2*n*k.Nt+jj,ii)  = sum(ktp1.Mw.*xa);
            Sol.rho(2*n*k.Nt+jj,ii)    = Zfac_struct.c.*Sol.M_tot(2*n*k.Nt+jj,ii);
            Sol.rho_vap(2*n*k.Nt+jj,ii)= Sol.rho(2*n*k.Nt+jj,ii)*Sol.wa(2*n*k.Nt+jj,ii);
            Sol.rho_dry(2*n*k.Nt+jj,ii)= Sol.rho(2*n*k.Nt+jj,ii)-Sol.rho_vap(2*n*k.Nt+jj,ii);
        end
        F_dry       = k.A_breathing.*sin(pi./k.T_breathing.*Sol.t(2*n*k.Nt+jj))*ones(size(k.zz));
        Sol.va(2*n*k.Nt+jj,:)   = F_dry./(Sol.rho_dry(2*n*k.Nt+jj,:).*A_a);
        Sol.F_vap(2*n*k.Nt+jj,:)= A_a.*Sol.va(2*n*k.Nt+jj,:).*Sol.rho_vap(2*n*k.Nt+jj,:);
        Sol.Fa(2*n*k.Nt+jj,:)   = A_a.*Sol.va(2*n*k.Nt+jj,:).*Sol.rho(2*n*k.Nt+jj,:);
        Sol.vart(2*n*k.Nt+jj,:) = -k.F_b./(k.rho_b.*A_art);
        Sol.vven(2*n*k.Nt+jj,:) = k.F_b./(k.rho_b.*A_ven);
    end
    
    % Hydraulic diameters
    D_a      = 4.*A_a./gamma_a;
    D_art    = 4.*A_art./gamma_art;
    D_ven    = 4.*A_ven./gamma_ven;
    % Resistivity coefficients
    for j = 2*n*k.Nt+1:2*n*k.Nt+k.Nt
        for ii= 1:k.Ns
            xa          = [Sol.x_a(j,ii); (1-Sol.x_a(j,ii))*0.79; (1-Sol.x_a(j,ii))*0.21];
            R_am        = CALC_R_a_m(Sol.Ta(j,ii),xa,Sol.Tw(j,ii),D_a(ii),Sol.va(j,ii));
            Sol.R_aw(j,ii)  = R_am(1,1);
            R               = CALC_R_ij(D_art(ii),D_ven(ii));
            Sol.R_wm(j,ii)  = R.m_it;
            Sol.R_art(j,ii) = R.it_art;
            Sol.R_ven(j,ii) = R.it_ven;
            Sol.R_mu(j,ii)  = R_am(2,2);
            Sol.R_qmu(j,ii)  = R_am(1,2);
        end
    end
    
    % Fluxes
    for j = 2*n*k.Nt+1:2*n*k.Nt+k.Nt
        for ii= 1:k.Ns
            xa          = [Sol.x_a(j,ii); (1-Sol.x_a(j,ii))*0.79; (1-Sol.x_a(j,ii))*0.21];
            Ma          = sum(xa.*ktp1.Mw);
            mu_a        = TP_ChemicalPotential_a(Sol.Ta(j,ii),k.pa,xa,1);
            h_id        = IG_H_a(Sol.Ta(j,ii),xa);
            h_res       = CB_Enthalpy_a(Sol.Ta(j,ii),k.pa,xa);
            h_i         = (h_id.dZ + h_res.dN)./ktp1.Mw;
            h_a         = (h_id.H + h_res.H)./Ma;
            mu_m        = TP_ChemicalPotential_w(Sol.Tw(j,ii),k.pa,1);
            Sol.h_a(j,ii) = h_a;
            Sol.h_wa(j,ii)= h_i(1);
            h_id        = IG_H_w(Sol.Tw(j,ii));
            h_res       = CB_Enthalpy_w(Sol.Tw(j,ii),k.pa);
            Sol.h_w(j,ii)= (h_id.H + h_res.H)./ktp2.Mw;
            X                  = [1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii); -(mu_m.mu(1)./Sol.Tw(j,ii)-mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii))];
            R_am        = CALC_R_a_m(Sol.Ta(j,ii),xa,Sol.Tw(j,ii),D_a(ii),Sol.va(j,ii));
            JJ          = R_am\X;
            Sol.X_T(j,ii)      = X(1);
            Sol.X_mu(j,ii)     = X(2);
            Sol.Jw(j,ii)       = JJ(2);
            Sol.Jq_a(j,ii)     = JJ(1);      % [J/m2s]
            Sol.Jq_w(j,ii)     = (1./Sol.Tm(j,ii)  -1./Sol.Tw(j,ii))./Sol.R_wm(j,ii) ;      % [J/m2s]
            Sol.Jq_art(j,ii)   = (1./Sol.Tart(j,ii)-1./Sol.Tm(j,ii))./Sol.R_art(j,ii);      % [J/m2s]
            Sol.Jq_ven(j,ii)   = (1./Sol.Tven(j,ii)-1./Sol.Tm(j,ii))./Sol.R_ven(j,ii);
            Sol.sigma(j,ii)    = (1./Sol.Tw(j,ii)  -1./Sol.Ta(j,ii)).*Sol.Jq_a(j,ii).*gamma_a(ii)+...
                (cos_m(ii)./Sol.Tm(j,ii)  -cos_w(ii)./Sol.Tw(j,ii)).*Sol.Jq_w(j,ii).*gamma_a(ii)+...
                (cos_art(ii)./Sol.Tart(j,ii)-cos_m(ii)./Sol.Tm(j,ii)).*Sol.Jq_art(j,ii).*gamma_art(ii)+...
                (cos_ven(ii)./Sol.Tven(j,ii)-cos_m(ii)./Sol.Tm(j,ii)).*Sol.Jq_ven(j,ii).*gamma_ven(ii)+...
                (-(cos_w(ii).*mu_m.mu(1)./Sol.Tw(j,ii)-cos_a(ii).*mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(cos_w(ii)./Sol.Tw(j,ii)-cos_a(ii)./Sol.Ta(j,ii))).*Sol.Jw(j,ii).*gamma_a(ii);
        end
        
    end

    %% Second half of the cycle

     M     = eye(5*k.Ns+k.Ns_art,5*k.Ns+k.Ns_art);
     M(1,1)= 0;
     M(k.Ns*1+1,k.Ns*1+1)= 0;
     M(k.Ns*2+1,k.Ns*2+1)= 0;
     M(k.Ns*5+1,k.Ns*5+1)= 0;
     M(k.Ns*4+1,k.Ns*4+1)= 0;
     M(k.Ns*4,k.Ns*4)= 0;
     Options = odeset('Mass',M);
     tt      = linspace(2.*n.*k.T_breathing+k.T_breathing,2.*n.*k.T_breathing+2*k.T_breathing,k.Nt);
     
     % Initial conditions
     T0a     = fliplr(Sol.Ta(end,:))./k.Scales.T;
     T0w     = fliplr(Sol.Tw(end,:))./k.Scales.T;
     T0m     = fliplr(Sol.Tm(end,:))./k.Scales.T;
     T0art   = fliplr(Sol.Tart_art(end,:))./k.Scales.T;
     T0ven   = fliplr(Sol.Tven(end,:))./k.Scales.T;
     w0a     = fliplr(Sol.wa(end,:));
     
     y0  = [ T0a T0w T0m T0ven w0a T0art];
     
     [t, x]  = ode15s(@DIFF_OnlyAir2, tt./k.Scales.t, y0,Options);
     disp("OUT == ")
    disp(t)
    disp(x)     
     % Time
     Sol.t(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt) = t.*k.Scales.t;
     
     % Tempeartures
     Sol.Ta(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,1:k.Ns)).*k.Scales.T;
     Sol.Tw(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,k.Ns+1:2*k.Ns)).*k.Scales.T;
     Sol.Tm(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,2*k.Ns+1:3*k.Ns)).*k.Scales.T;
     Sol.Tven(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)= fliplr(x(:,3*k.Ns+1:4*k.Ns)).*k.Scales.T;
     Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)= fliplr(x(:,4*k.Ns+1:5*k.Ns));
     Sol.x_a(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:).*k.M_dry./(ktp1.Mw(1).*(1-Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:))+Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:).*k.M_dry);
%     Sol.x_a(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:) = (0.21.*Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:).*ktp1.Mw(3)./ktp1.Mw(1)+0.79.*Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:).*ktp1.Mw(2)./ktp1.Mw(1))./(1-Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)+0.21.*Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:).*ktp1.Mw(3)./ktp1.Mw(1)+0.89.*Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:).*ktp1.Mw(2)./ktp1.Mw(1));
     Sol.Tart_art(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)= fliplr(x(:,5*k.Ns+1:5*k.Ns+k.Ns_art)).*k.Scales.T;
     Sol.Tart(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)= pchip(k.zz_art,Sol.Tart_art(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:),Sol.x);
     for jj = k.Nt+1:2*k.Nt
        for ii= 1:k.Ns
            xa          = [Sol.x_a(jj,ii); (1-Sol.x_a(jj,ii))*0.79; (1-Sol.x_a(jj,ii))*0.21];
            [Zfac_struct,~] = TP_Zfac_a(Sol.Ta(jj,ii),k.pa,xa,2);
            Sol.M_tot(2*n*k.Nt+jj,ii)  = sum(ktp1.Mw.*xa);
            Sol.rho(2*n*k.Nt+jj,ii)    = Zfac_struct.c.*Sol.M_tot(2*n*k.Nt+jj,ii);
            Sol.rho_vap(2*n*k.Nt+jj,ii)= Sol.rho(2*n*k.Nt+jj,ii)*Sol.wa(2*n*k.Nt+jj,ii);
            Sol.rho_dry(2*n*k.Nt+jj,ii)= Sol.rho(2*n*k.Nt+jj,ii)-Sol.rho_vap(2*n*k.Nt+jj,ii);
        end
        F_dry       = k.A_breathing.*sin(pi./k.T_breathing.*Sol.t(2*n*k.Nt+jj))*ones(size(k.zz));
        Sol.va(2*n*k.Nt+jj,:)   = F_dry./(Sol.rho_dry(2*n*k.Nt+jj,:).*A_a);
        Sol.F_vap(2*n*k.Nt+jj,:)= A_a.*Sol.va(2*n*k.Nt+jj,:).*Sol.rho_vap(2*n*k.Nt+jj,:);
        Sol.Fa(2*n*k.Nt+jj,:)   = A_a.*Sol.va(2*n*k.Nt+jj,:).*Sol.rho(2*n*k.Nt+jj,:);
        Sol.vart(2*n*k.Nt+jj,:) = -k.F_b./(k.rho_b.*A_art);
        Sol.vven(2*n*k.Nt+jj,:) = k.F_b./(k.rho_b.*A_ven);
    end
   
     
     % Hydraulic diameters
     D_a      = 4.*A_a./gamma_a;
     D_art    = 4.*A_art./gamma_art;
     D_ven    = 4.*A_ven./gamma_ven;
    
    % Resistivity coefficients
        for j = 2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt
            for ii= 1:k.Ns
                xa          = [Sol.x_a(j,ii); (1-Sol.x_a(j,ii))*0.79; (1-Sol.x_a(j,ii))*0.21];
                R_am        = CALC_R_a_m(Sol.Ta(j,ii),xa,Sol.Tw(j,ii),D_a(ii),Sol.va(j,ii));
                Sol.R_aw(j,ii)  = R_am(1,1);
                R               = CALC_R_ij(D_art(ii),D_ven(ii));
                Sol.R_wm(j,ii)  = R.m_it;
                Sol.R_art(j,ii) = R.it_art;
                Sol.R_ven(j,ii) = R.it_ven;
                Sol.R_mu(j,ii)  = R_am(2,2);
            Sol.R_qmu(j,ii)  = R_am(1,2);
            end
        end
    % Fluxes
    for j = 2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt
        for ii= 1:k.Ns
            xa          = [Sol.x_a(j,ii); (1-Sol.x_a(j,ii))*0.79; (1-Sol.x_a(j,ii))*0.21];
            Ma          = sum(xa.*ktp1.Mw);
            mu_a        = TP_ChemicalPotential_a(Sol.Ta(j,ii),k.pa,xa,1);
            h_id        = IG_H_a(Sol.Ta(j,ii),xa);
            h_res       = CB_Enthalpy_a(Sol.Ta(j,ii),k.pa,xa);
            h_i         = (h_id.dZ + h_res.dN)./ktp1.Mw;
            h_a         = (h_id.H + h_res.H)./Ma;
            mu_m        = TP_ChemicalPotential_w(Sol.Tw(j,ii),k.pa,1);
            Sol.h_a(j,ii) = h_a;
            Sol.h_wa(j,ii)= h_i(1);
            h_id        = IG_H_w(Sol.Tw(j,ii));
            h_res       = CB_Enthalpy_w(Sol.Tw(j,ii),k.pa);
            Sol.h_w(j,ii)= (h_id.H + h_res.H)./ktp2.Mw;
            X                  = [1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii); -(mu_m.mu(1)./Sol.Tw(j,ii)-mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(1./Sol.Tw(j,ii)-1./Sol.Ta(j,ii))];
            R_am        = CALC_R_a_m(Sol.Ta(j,ii),xa,Sol.Tw(j,ii),D_a(ii),Sol.va(j,ii));
            JJ          = R_am\X;
            Sol.X_T(j,ii)      = X(1);
            Sol.X_mu(j,ii)     = X(2);
            Sol.Jw(j,ii)       = JJ(2);
            Sol.Jq_a(j,ii)     = JJ(1);      % [J/m2s]
            Sol.Jq_w(j,ii)     = (1./Sol.Tm(j,ii)  -1./Sol.Tw(j,ii))./Sol.R_wm(j,ii) ;      % [J/m2s]
            Sol.Jq_art(j,ii)   = (1./Sol.Tart(j,ii)-1./Sol.Tm(j,ii))./Sol.R_art(j,ii);      % [J/m2s]
            Sol.Jq_ven(j,ii)   = (1./Sol.Tven(j,ii)-1./Sol.Tm(j,ii))./Sol.R_ven(j,ii);
            Sol.sigma(j,ii)    = (1./Sol.Tw(j,ii)  -1./Sol.Ta(j,ii)).*Sol.Jq_a(j,ii).*gamma_a(ii)+...
                (cos_m(ii)./Sol.Tm(j,ii)  -cos_w(ii)./Sol.Tw(j,ii)).*Sol.Jq_w(j,ii).*gamma_a(ii)+...
                (cos_art(ii)./Sol.Tart(j,ii)-cos_m(ii)./Sol.Tm(j,ii)).*Sol.Jq_art(j,ii).*gamma_art(ii)+...
                (cos_ven(ii)./Sol.Tven(j,ii)-cos_m(ii)./Sol.Tm(j,ii)).*Sol.Jq_ven(j,ii).*gamma_ven(ii)+...
                (-(cos_w(ii).*mu_m.mu(1)./Sol.Tw(j,ii)-cos_a(ii).*mu_a.mu(1)./Sol.Ta(j,ii))+h_i(1).*(cos_w(ii)./Sol.Tw(j,ii)-cos_a(ii)./Sol.Ta(j,ii))).*Sol.Jw(j,ii).*gamma_a(ii);
        end
    end
    
   if i == 2 && k.initCond == 0 
       resTexp0  = abs(Sol.Ta(end - k.Nt/2,1)         - Sol.Ta(end - k.Nt/2 - 2*k.Nt,1)) ;
       resT0     = norm(Sol.Ta(end-2*k.Nt+1:end,:)    - Sol.Ta(end-4*k.Nt+1:end-2*k.Nt,:)) ;
       resSigma0 = norm(Sol.sigma(end-2*k.Nt+1:end,:) - Sol.sigma(end-4*k.Nt+1:end-2*k.Nt,:)) ;
       
   elseif i == 1 && k.initCond ~= 0 
       Nose0 = load(k.initFileName);
       resTexp0  = abs(Nose0.Nose.Ta(4*k.Nt-k.Nt/2,1)       - Nose0.Nose.Ta(2*k.Nt - k.Nt/2 ,1)) ;
       resT0     = norm(Nose0.Nose.Ta(2*k.Nt+1:4*k.Nt,:)    - Nose0.Nose.Ta(1:2*k.Nt,:)) ;
       resSigma0 = norm(Nose0.Nose.sigma(2*k.Nt+1:4*k.Nt,:) - Nose0.Nose.sigma(1:2*k.Nt,:)) ;
       
   end
  
       
   
   
   if (i > 2 && k.initCond == 0) || (i > 1 && k.initCond ~= 0) 
       resTexp   = abs(Sol.Ta(end - k.Nt/2,1)         - Sol.Ta(end - k.Nt/2 - 2*k.Nt,1)) ;
       resT      = norm(Sol.Ta(end-2*k.Nt+1:end,:)    - Sol.Ta(end-4*k.Nt+1:end-2*k.Nt,:)) ;       
       resSigma  = norm(Sol.sigma(end-2*k.Nt+1:end,:) - Sol.sigma(end-4*k.Nt+1:end-2*k.Nt,:)) ;

       disp("resTexp0 == ")
       disp(resTexp0)
       disp("resT0 ==")
       disp(resT0)
       disp("resSigma0 ==")
       disp(resSigma0)
       
       disp("== resTexp == ")
       disp(resTexp)
       disp("== resT ==")
       disp(resT)
       disp("== resSigma ==")
       disp(resSigma)       
       
       disp("== resTexp/resTexp0 == ")
       disp(resTexp/resTexp0)
       disp("== resT/resT0 ==")
       disp(resT/resT0)
       disp("== resSigma/resSigma0 ==")
       disp(resSigma/resSigma0)
       
       for j = 1:k.Ns
           totSigma(j) = trapz(Sol.t(end-k.Nt*2+1:end),    Sol.sigma(end-k.Nt*2+1:end,j)); 
           if totSigma(j) < 0
               disp("Negative local entropy production. ")
               disp("j")
               disp(totSigma(j)) 
           end
       end
       
       
       if resTexp/resTexp0 < 1e-5
           disp("Expelled temperature is converged. 1e-5 ")
           disp("N ==  ")
           disp(i)
       end

       if resT/resT0 < 1e-6
           disp("Norm Air temperature is converged. 1e-6")
           %disp("N ==  ")
           %disp(i)
           %break;
       end

       if resSigma/resSigma0 < 1e-6
           disp("Norm Sigma (entropy) is converged. 1e-6")
           %disp("N ==  ")
           %disp(i)
           %break;
       end
       
       if resTexp/resTexp0 < 1e-5 && resT/resT0 < 1e-6 && resSigma/resSigma0 < 1e-6
           disp("All converged.! ")
           
           if k.convCond == 1
               disp(" FINISH iteration")
               break;
           end
           
       end
   
   end
    
    
    
    
    
        % Initial conditions
        T0a     = (Sol.Ta(end,:))./k.Scales.T;
        T0w     = (Sol.Tw(end,:))./k.Scales.T;
        T0m     = (Sol.Tm(end,:))./k.Scales.T;
        T0art   = (Sol.Tart_art(end,:))./k.Scales.T;
        T0ven   = (Sol.Tven(end,:))./k.Scales.T;
        w0a     = (Sol.wa(end,:));
       y0  = [ T0a T0w T0m  T0ven w0a T0art];
       
       Nose = Sol;
       save(k.fileName, 'Nose');
       %save RD_Ns36_36_veindssrev_1000c_m30.mat Sol
       %save RD_Ns24_24_veindssrev_1000c_m30.mat Sol
       %save SA_shunt_Ns24_24_veindssrev_1000c_m30.mat Sol
end
