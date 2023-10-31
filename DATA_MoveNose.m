function sol = DATA_MoveNose(fileNameA, fileNameB, fileNameC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% === input === 
% fileNameA: str, file directory, move into this file from the index 1. 
% fileNameB: str, file directory, move from this file to fileName1.
% fileNameC: str, file name for save. 
%
% === output ===
% sol.(variable below)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noseA = load(fileNameA,'Nose'); 
noseB = load(fileNameB,'Nose'); 

 
noseC.Nose.t        = cat(2, noseA.Nose.t, noseB.Nose.t + noseA.Nose.t(end));
if length(noseA.Nose.x) ~= length(noseB.Nose.x)
   disp("ERROR! x-scales are different! ") 
end
noseC.Nose.x        = noseA.Nose.x;
noseC.Nose.Ta       = cat(1, noseA.Nose.Ta,       noseB.Nose.Ta);
noseC.Nose.Tw       = cat(1, noseA.Nose.Tw,       noseB.Nose.Tw);
noseC.Nose.Tm       = cat(1, noseA.Nose.Tm,       noseB.Nose.Tm);
noseC.Nose.Tven     = cat(1, noseA.Nose.Tven,     noseB.Nose.Tven);
noseC.Nose.wa       = cat(1, noseA.Nose.wa,       noseB.Nose.wa);
noseC.Nose.x_a      = cat(1, noseA.Nose.x_a,      noseB.Nose.x_a);
noseC.Nose.Tart_art = cat(1, noseA.Nose.Tart_art, noseB.Nose.Tart_art);
noseC.Nose.Tart     = cat(1, noseA.Nose.Tart,     noseB.Nose.Tart);
noseC.Nose.M_tot    = cat(1, noseA.Nose.M_tot,    noseB.Nose.M_tot);
noseC.Nose.rho      = cat(1, noseA.Nose.rho,      noseB.Nose.rho);
noseC.Nose.rho_vap  = cat(1, noseA.Nose.rho_vap,  noseB.Nose.rho_vap);
noseC.Nose.rho_dry  = cat(1, noseA.Nose.rho_dry,  noseB.Nose.rho_dry);
noseC.Nose.va       = cat(1, noseA.Nose.va,       noseB.Nose.va);
noseC.Nose.F_vap    = cat(1, noseA.Nose.F_vap,    noseB.Nose.F_vap);
noseC.Nose.Fa       = cat(1, noseA.Nose.Fa,       noseB.Nose.Fa);
noseC.Nose.vart     = cat(1, noseA.Nose.vart,     noseB.Nose.vart);
noseC.Nose.vven     = cat(1, noseA.Nose.vven,     noseB.Nose.vven);
noseC.Nose.R_aw     = cat(1, noseA.Nose.R_aw,     noseB.Nose.R_aw);
noseC.Nose.R_wm     = cat(1, noseA.Nose.R_wm,     noseB.Nose.R_wm);
noseC.Nose.R_art    = cat(1, noseA.Nose.R_art,    noseB.Nose.R_art);
noseC.Nose.R_ven    = cat(1, noseA.Nose.R_ven,    noseB.Nose.R_ven);
noseC.Nose.R_mu     = cat(1, noseA.Nose.R_mu,     noseB.Nose.R_mu);
noseC.Nose.R_qmu    = cat(1, noseA.Nose.R_qmu,    noseB.Nose.R_qmu );
noseC.Nose.h_a      = cat(1, noseA.Nose.h_a,      noseB.Nose.h_a);
noseC.Nose.h_wa     = cat(1, noseA.Nose.h_wa,     noseB.Nose.h_wa);
noseC.Nose.h_w      = cat(1, noseA.Nose.h_w,      noseB.Nose.h_w);
noseC.Nose.X_T      = cat(1, noseA.Nose.X_T,      noseB.Nose.X_T);
noseC.Nose.X_mu     = cat(1, noseA.Nose.X_mu,     noseB.Nose.X_mu);

noseC.Nose.Jw       = cat(1, noseA.Nose.Jw,       noseB.Nose.Jw);
noseC.Nose.Jq_a     = cat(1, noseA.Nose.Jq_a,     noseB.Nose.Jq_a);
noseC.Nose.Jq_w     = cat(1, noseA.Nose.Jq_w,     noseB.Nose.Jq_w);
noseC.Nose.Jq_art   = cat(1, noseA.Nose.Jq_art,   noseB.Nose.Jq_art);
noseC.Nose.Jq_ven   = cat(1, noseA.Nose.Jq_ven,   noseB.Nose.Jq_ven);
noseC.Nose.sigma    = cat(1, noseA.Nose.sigma,    noseB.Nose.sigma);

%save noseC 

%save fileName.mat noseC.Nose % Save the results as a matrix
Nose = noseC.Nose;

save(fileNameC, 'Nose');
%fileName = [fileNameC, '.mat'];

%save fileName sol

