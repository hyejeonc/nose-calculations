function R = CALC_R_ij(D_art,D_ven)

global k

delta_vessel = 3e-5; % 0.2mm for reindeer, but 200 um for seals 
h_art   = 3.66*k.k_b./D_art;              
h_ven   = 3.66*k.k_b./D_ven; 

R.m_it  = k.factor_Rij*(k.d_mucus/(k.k_w*k.T_body^2)*100);%+ k.delta_tissue/(k.k_b*k.T_body^2)*10);
R.it_art= (delta_vessel/(k.k_b*k.T_body^2)+ 1/(h_art*k.T_body^2))*10;
R.it_ven= (delta_vessel/(k.k_b*k.T_body^2)+ 1/(h_ven*k.T_body^2))*10;

