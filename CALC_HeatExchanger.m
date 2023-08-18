function Sol = CALC_HeatExchanger

global k ktp1 ktp2

k = ScalesParameters(k);

wa      = k.wa;
wb      = k.wb;


% Initial conditions
% T0a     = linspace(k.Ta,k.T_body,k.Ns)./k.Scales.T;
% T0w     = linspace(k.T_body,k.T_body,k.Ns)./k.Scales.T;
% T0m     = linspace(k.T_body,k.T_body,k.Ns)./k.Scales.T;
% T0art   = linspace(k.T_body,k.T_body,k.Ns_art)./k.Scales.T;
% T0ven   = linspace(k.T_body,k.T_body,k.Ns)./k.Scales.T;
% w0a     = linspace(wa,wb,k.Ns);
% y0  = [ T0a T0w T0m  T0ven w0a T0art];

load('Ref_0.mat');
T0a     = Ref_Case.Ta(end-2*k.Nt+1,:)./k.Scales.T;
T0w     = Ref_Case.Tw(end-2*k.Nt+1,:)./k.Scales.T;
T0m     = Ref_Case.Tm(end-2*k.Nt+1,:)./k.Scales.T;
T0art   = pchip(k.zz,Ref_Case.Tart(end-2*k.Nt+1,:),k.zz_art)./k.Scales.T;
T0ven   = Ref_Case.Tven(end-2*k.Nt+1,:)./k.Scales.T;
w0a     = Ref_Case.wa(end-2*k.Nt+1,:);
y0  = [ T0a T0w T0m  T0ven w0a T0art];


for i = 1:k.N_cycle
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
    
    [t, x] = ode15s(@DIFF_HeatExchanger,tt./k.Scales.t,y0,Options);
    
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
    
    gamma_a  = k.P_a_mean*ones(size(k.zz));
    gamma_art= k.P_art_mean*ones(size(k.zz));
    gamma_ven= k.P_ven_mean*ones(size(k.zz));
    
    A_a      = k.A_a_mean*ones(size(k.zz));    % [m2]
    A_art    = k.A_art_mean*ones(size(k.zz));     % [m2]
    A_ven    = k.A_ven_mean*ones(size(k.zz));     % [m2]
    
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
            Sol.X_T(j,ii)  = X(1);
            Sol.X_mu(j,ii) = X(2);
            Sol.Jw(j,ii)       = JJ(2);
            Sol.Jq_a(j,ii)     = JJ(1);      % [J/m2s]
            Sol.Jq_w(j,ii)     = (1./Sol.Tm(j,ii)  -1./Sol.Tw(j,ii))./Sol.R_wm(j,ii) ;      % [J/m2s]
            Sol.Jq_art(j,ii)   = (1./Sol.Tart(j,ii)-1./Sol.Tm(j,ii))./Sol.R_art(j,ii);      % [J/m2s]
            Sol.Jq_ven(j,ii)   = (1./Sol.Tven(j,ii)-1./Sol.Tm(j,ii))./Sol.R_ven(j,ii);
            Sol.sigma(j,ii)    = (1./Sol.Tw(j,ii)  -1./Sol.Ta(j,ii)).*Sol.Jq_a(j,ii).*gamma_a(ii)+...
                (1./Sol.Tm(j,ii)  -1./Sol.Tw(j,ii)).*Sol.Jq_w(j,ii).*gamma_a(ii)+...
                (1./Sol.Tart(j,ii)-1./Sol.Tm(j,ii)).*Sol.Jq_art(j,ii).*gamma_art(ii)+...
                (1./Sol.Tven(j,ii)-1./Sol.Tm(j,ii)).*Sol.Jq_ven(j,ii).*gamma_ven(ii)+...
                X(2).*Sol.Jw(j,ii).*gamma_a(ii);
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
    
    [t, x]  = ode15s(@DIFF_HeatExchanger2, tt./k.Scales.t, y0,Options);
    
    % Time
    Sol.t(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt) = t.*k.Scales.t;
    
    % Tempeartures
    Sol.Ta(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,1:k.Ns)).*k.Scales.T;
    Sol.Tw(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,k.Ns+1:2*k.Ns)).*k.Scales.T;
    Sol.Tm(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,2*k.Ns+1:3*k.Ns)).*k.Scales.T;
    Sol.Tven(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)= fliplr(x(:,3*k.Ns+1:4*k.Ns)).*k.Scales.T;
    Sol.wa(2*n*k.Nt+k.Nt+1:2*n*k.Nt+2*k.Nt,:)  = fliplr(x(:,4*k.Ns+1:5*k.Ns));
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
    %
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
            Sol.X_T(j,ii)  = X(1);
            Sol.X_mu(j,ii) = X(2);
            Sol.Jw(j,ii)       = JJ(2);
            Sol.Jq_a(j,ii)     = JJ(1);      % [J/m2s]
            Sol.Jq_w(j,ii)     = (1./Sol.Tm(j,ii)  -1./Sol.Tw(j,ii))./Sol.R_wm(j,ii) ;      % [J/m2s]
            Sol.Jq_art(j,ii)   = (1./Sol.Tart(j,ii)-1./Sol.Tm(j,ii))./Sol.R_art(j,ii);      % [J/m2s]
            Sol.Jq_ven(j,ii)   = (1./Sol.Tven(j,ii)-1./Sol.Tm(j,ii))./Sol.R_ven(j,ii);
            Sol.sigma(j,ii)    = (1./Sol.Tw(j,ii)  -1./Sol.Ta(j,ii)).*Sol.Jq_a(j,ii).*gamma_a(ii)+...
                (1./Sol.Tm(j,ii)  -1./Sol.Tw(j,ii)).*Sol.Jq_w(j,ii).*gamma_a(ii)+...
                (1./Sol.Tart(j,ii)-1./Sol.Tm(j,ii)).*Sol.Jq_art(j,ii).*gamma_art(ii)+...
                (1./Sol.Tven(j,ii)-1./Sol.Tm(j,ii)).*Sol.Jq_ven(j,ii).*gamma_ven(ii)+...
                X(2).*Sol.Jw(j,ii).*gamma_a(ii);
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
end