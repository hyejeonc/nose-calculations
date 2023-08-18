function R = CALC_R_i

global k

R.a  = 1/(k.k_a*k.T_ref^2);
R.m  = 1/(k.k_w*k.T_ref^2);
R.it = 1/(k.k_b*k.T_ref^2);
R.art= 1/(k.k_b*k.T_body^2);
R.ven= 1/(k.k_b*k.T_body^2);