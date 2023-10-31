function result = DATA_Txtdata(n)

close all

% This returns position, perimeter and area data. 
% xi, yiPeri and yiArea all have n x 1 dimension. 
%=============================================================
% input
% n = integer
% result.xArcHuntley    = x_nscale     %% x_nscaled ; 
% result.periArcHuntley = peri_nscale  %% perimeter_nscaled ;
% result.areaArcHuntley = area_nscale  %% area_nscaled; 
%
% From Huntley et al, 1984
%
% ============================================================
% Arctic = Bearded seal (Eb)
% Med = Mediterranean seal (Mm)
% 1 paper from Mason
% 2 data for Eirik, excluding direct channel to olfactory 
% 3 data for HJ, including direct channel to olfactory 

% M = readtable('sealdata.csv');%,'NumHeaderLines',1);

% result.xArcHuntley    = x_nscale; 
% result.periArcHuntley = peri_nscale;
% result.areaArcHuntley = area_nscale; 


% result.xArc = xiArc;
% result.xArcPosi = xArcPosi; 
% result.xArcPeri = xArcPeri ;
% result.xArcAreaWshunt = xArcAreaWshunt;
% result.xArcAreaWoshunt = xArcAreaWoshunt;
% 
% result.xSub = xiSub; 
% result.xSubPosi = xSubPosi; 
% result.xSubPeri = xSubPeri ;
% result.xSubAreaWshunt = xSubAreaWshunt;
% result.xSubAreaWoshunt = xSubAreaWoshunt;

resMason = DATA_Csvdata(36, 0.5);
%disp(resMason.xArc)
disp("resMason.xArcPosi")
disp(max(resMason.xArcPosi) - min(resMason.xArcPosi))
%disp("resMason.xArcPeri")
%disp(max(resMason.xArcPeri) - min(resMason.xArcPeri))
%disp("resMason.xArcAreaWshunt")
%disp(min(resMason.xArcAreaWshunt) - max(resMason.xArcAreaWshunt))
disp("L == ")
disp(resMason.xArc)
%disp(resMason.xArcPosi(end)*resMason.xArc(2) - resMason.xArcPosi(end)*resMason.xArc(1));
%disp(max(resMason.xArcPosi) - max(resMason.xArcPosi));


M = readtable('elephant_seal_data_huntley1984.csv');%,'NumHeaderLines',1);

x    = table2array(M(1:end,1)) ; % perimeter 
peri = table2array(M(1:end,2)) ; % perimeter 
area = table2array(M(1:end,3)) ; % perimeter 

%disp("x")
%disp(max(x) - min(x))
%disp("arctic_peri")
%disp(peri)
%disp("aricicarea")
%disp(area)

x_nscale    = linspace(resMason.xArc(2)*max(x), resMason.xArc(1)*max(x), n);
peri_nscale = pchip(x, peri, x_nscale);
area_nscale = pchip(x, area, x_nscale);


%disp("x_scale")
%disp(x_nscale)
%disp("peri_scale")
%disp(peri_nscale)
%disp("area_scale")
%disp(area_nscale)
% xArcPosi        = linspace(xiArc(1)*max(A(:,2)), xiArc(2)*max(A(:,2)), n); 
% 
% 
% xSubPosi        = linspace(xiSub(1)*max(S(:,7)), xiSub(2)*max(S(:,7)), n); 
% 
% xSubPeri        = pchip(S(:,7), S(:,8),  xSubPosi);
% xSubAreaWoshunt = pchip(S(:,7), S(:,10), xSubPosi);
% xSubAreaWshunt  = pchip(S(:,7), S(:,11), xSubPosi);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %[xiArc,yiArc] = polyxpoly(x/max(x), A(:,5)./A(:,6), A(:,2)/max(A(:,2)), ratio.*ones(size(A(:,2))) );
% 
% 
% 
% %%%%%%%%%%%%%%%
% %%% from DATA_Csvdata.m 
% 
% A(:,2) = table2array(M(1:end-8,2)) ; % arctic seal position (*)
% A(:,3) = table2array(M(1:end-8,3)) ; % arctic seal perimeter (*)
% A(:,4) = table2array(M(1:end-8,4)) ; % arctic seal area + only bone
% A(:,5) = table2array(M(1:end-8,5)) ; % arctic seal area + w/o shunt *** 
% A(:,6) = table2array(M(1:end-8,6)) ; % arctic seal area + w/ shunt *** (*)
% 
% S(:,7) = table2array(M(:,7)) ; % subtropical seal position (*)
% S(:,8) = table2array(M(:,8)) ; % subtropical seal perimeter (*)
% S(:,9) = table2array(M(:,9)) ; % subtropical seal area + only bone
% S(:,10) = table2array(M(:,10)) ; % subtropical seal area + w/o shunt ***
% S(:,11) = table2array(M(:,11)) ; % subtropical seal area + w/ shunt *** (*)
% 
% %ratio = 0.5;
% %n = 36;
% % I wanna see proportion of area with shunt divided by area without shunt 
% % figure(1)
% % plot(A(:,2)/max(A(:,2)), A(:,5)./A(:,6), 'DisplayName','Arctic' )  ; hold on
% % plot(A(:,7)/max(A(:,7)), A(:,10)./A(:,11), 'DisplayName','Subtropical' )  ; hold on
% % plot(A(:,2)/max(A(:,2)), ratio.*ones(size(A(:,2)/max(A(:,2)))), 'DisplayName', 'y=0.5' )  ;
% % xlabel('z/L [-]','FontSize',14);
% % ylabel('A_{MT} / A_{MT,shunt} [-]','FontSize',14);
% % legend('FontSize',14, 'Location','south'); hold off 
% % saveas(figure(1), [pwd '/fig/areaperi_divide.pdf']);
% % saveas(figure(1), [pwd '/fig/areaperi_divide.eps'], 'epsc');
% 
% % disp(A(:,2)/max(A(:,2))) ;
% % disp(A(:,5)./A(:,6)) ;
% % disp(A(:,2)/max(A(:,2))) ;
% % disp(ratio.*ones(size(A(:,2)))) ;
% % [xi,yi] = polyxpoly(A(:,2)/max(A(:,2)), A(:,5)./A(:,6),A(:,7)/max(A(:,7)), 0.5.*ones(size(A(:,7)/max(A(:,7)))));
% 
% 
% [xiArc,yiArc] = polyxpoly(A(:,2)/max(A(:,2)), A(:,5)./A(:,6), A(:,2)/max(A(:,2)), ratio.*ones(size(A(:,2))) );
% % find cross section point [xi(1), yi(1)], [xi(2), yi(2)] for arctic seal. 
% 
% % disp("xiArc")
% % disp(xiArc)
% % disp("yiArc")
% % disp(yiArc)
% 
% 
% 
% 
% %% Find the relation of dosal meatus area and specie  
% %figure(2) 
% %plot(A(:,2)/max(A(:,2)), A(:,6)-A(:,5), 'DisplayName','Arctic' )  ; hold on
% %plot(S(:,7)/max(S(:,7)), S(:,11)-S(:,10), 'DisplayName','Subtropical' )  ;
% %xlabel('z/L [-]','FontSize',14);
% %ylabel('A_{D} [mm^2]','FontSize',14);
% %legend('FontSize',14); hold off
% % saveas(figure(2),'./fig/areaperi_subtract.pdf')
% % saveas(figure(2),'./fig/areaperi_subtract.eps', 'epsc')
% 
% 
% %figure(3)
% %plot(A(:,2)/max(A(:,2)), A(:,6), 'DisplayName','Arctic' )  ; hold on
% %plot(S(:,7)/max(S(:,7)), S(:,11), 'DisplayName','Subtropical' )  ; 
% %%%plot(A(:,7)/A(end,7), A(:,11), 'DisplayName','Subtropical' )  ; 
% %scatter(xArcPosi/max(A(:,2)), xArcAreaWshunt, 'DisplayName','Arctic' )  ; 
% %scatter(xSubPosi/max(S(:,7)), xSubAreaWshunt, 'DisplayName','Subtropical' )  ; 
% % disp("what is A(end,2) ?? ")
% % disp(S(end,7))
% % disp(S(:,7)./S(end,7))
% % disp(S(:,11))
% %xlabel('z/L [-]','FontSize',14);
% %ylabel('A [mm^2]','FontSize',14);
% 
% %legend('Location', 'Best','FontSize',14); 
% %saveas(figure(3),'./fig/area_subtract_new.pdf')
% %saveas(figure(3),'./fig/area_subtract_new.eps', 'epsc')
% 
% 
% %figure(4)
% %plot(A(:,2)/max(A(:,2)), A(:,3), 'DisplayName','Arctic' )  ; hold on
% %plot(S(:,7)/max(S(:,7)), S(:,8), 'DisplayName','Subtropical' )  ; 
% %scatter(xArcPosi/max(A(:,2)), xArcPeri, 'DisplayName','Arctic' )  ; 
% %scatter(xSubPosi/max(S(:,7)), xSubPeri, 'DisplayName','Subtropical' )  ; 
% %xlabel('z/L [-]','FontSize',14);
% %ylabel('\gamma [mm]','FontSize',14);
% 
% %legend('Location', 'northwest','FontSize',14);
% %saveas(figure(4), [pwd '/fig/peri_subtract_new.pdf'])
% %saveas(figure(4),'./fig/peri_subtract_new.eps', 'epsc')
% 
% result.xArc = xiArc;
% result.xArcPosi = xArcPosi; 
% result.xArcPeri = xArcPeri ;
% result.xArcAreaWshunt = xArcAreaWshunt;
% result.xArcAreaWoshunt = xArcAreaWoshunt;
% 
% result.xSub = xiSub; 
% result.xSubPosi = xSubPosi; 
% result.xSubPeri = xSubPeri ;
% result.xSubAreaWshunt = xSubAreaWshunt;
% result.xSubAreaWoshunt = xSubAreaWoshunt;
% 
% %result = [xArc; xArcPeri; xArcAreaWshunt; xSub; xSubPeri; xSubAreaWshunt];
% 
% 
% xArcPosi        = linspace(xiArc(1)*max(A(:,2)), xiArc(2)*max(A(:,2)), n); 
% xArcPeri        = pchip(A(:,2), A(:,3), xArcPosi);
% xArcAreaWoshunt = pchip(A(:,2), A(:,5), xArcPosi);
% xArcAreaWshunt  = pchip(A(:,2), A(:,6), xArcPosi);
% 
% 
% 
% 
% 
% xSubPosi        = linspace(xiSub(1)*max(S(:,7)), xiSub(2)*max(S(:,7)), n); 
% 
% x_nscale       = linspace(xiArc(1)*max(A(:,2)), xiArc(2)*max(A(:,2)), n); 

result.xArcHuntley    = x_nscale ; 
result.periArcHuntley = peri_nscale ;
result.areaArcHuntley = area_nscale; 
%disp(x_nscale)
%disp(peri_nscale)
%disp(area_nscale)
end

