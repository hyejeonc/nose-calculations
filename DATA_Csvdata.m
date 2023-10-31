function result = DATA_Csvdata(n, ratio)
close all

% This returns position, perimeter and area data. 
% xi, yiPeri and yiArea all have n x 1 dimension. 
%=============================================================
% r0: find ratio by col K and N from excel(sealdata.csv) - area of air and area of air and shunt
% r1: find ratio by col K-H and col N-H 
%=============================================================
% INPUT
% n = integer
% ratio = float
%
%
% OUTPUT
% result.xArc = array, x intersection position for arctic seal ( 0 < result.xArc < 1). If the number of intersection is more than one, it is array. 
% result.xArcPosi = xArcPosi; 
% result.xArcPeri = xArcPeri;
% result.xArcAreaWshunt = xArcAreaWshunt;
% result.xArcAreaWoshunt = xArcAreaWoshunt;
% 
% result.xSub = xiSub; 
% result.xSubPosi = xSubPosi; 
% result.xSubPeri = xSubPeri ;
% result.xSubAreaWshunt = xSubAreaWshunt;
% result.xSubAreaWoshunt = xSubAreaWoshunt;
% ============================================================
% Arctic = Bearded seal (Eb)
% Med = Mediterranean seal (Mm)
% 1 paper from Mason
% 2 data for Eirik, excluding direct channel to olfactory 
% 3 data for HJ, including direct channel to olfactory 

M = readtable('sealdata.csv');%,'NumHeaderLines',1);

%% !! 5 - 4 !!!!! 
A(:,2) = table2array(M(1:end-8,2)) ; % arctic seal position (*)             col D
A(:,3) = table2array(M(1:end-8,3)) ; % arctic seal perimeter (*)            col F
A(:,4) = table2array(M(1:end-8,4)) ; % arctic seal area + only bone         col H
A(:,5) = table2array(M(1:end-8,5) ) ; % arctic seal area + w/o shunt ***    col K
A(:,6) = table2array(M(1:end-8,6)) ; % arctic seal area + w/ shunt *** (*)  col N

S(:,7) = table2array(M(:,7)) ; % subtropical seal position (*)              col D
S(:,8) = table2array(M(:,8)) ; % subtropical Arcseal perimeter (*)             col F
S(:,9) = table2array(M(:,9)) ; % subtropical seal area + only bone          col H        
S(:,10) = table2array(M(:,10)) ; % subtropical seal area + w/o shunt ***    col K
S(:,11) = table2array(M(:,11)) ; % subtropical seal area + w/ shunt *** (*) col N

%ratio = 0.5;
%n = 36;
% I wanna see proportion of area with shunt divided by area without shunt 
% figure(1)
% plot(A(:,2)/max(A(:,2)), A(:,5)./A(:,6), 'DisplayName','Arctic' )  ; hold on
% plot(A(:,7)/max(A(:,7)), A(:,10)./A(:,11), 'DisplayName','Subtropical' )  ; hold on
% plot(A(:,2)/max(A(:,2)), ratio.*ones(size(A(:,2)/max(A(:,2)))), 'DisplayName', 'y=0.5' )  ;
% xlabel('z/L [-]','FontSize',14);
% ylabel('A_{MT} / A_{MT,shunt} [-]','FontSize',14);
% legend('FontSize',14, 'Location','south'); hold off 
% saveas(figure(1), [pwd '/fig/areaperi_divide.pdf']);
% saveas(figure(1), [pwd '/fig/areaperi_divide.eps'], 'epsc');

% disp(A(:,2)/max(A(:,2))) ;
% disp(A(:,5)./A(:,6)) ;
% disp(A(:,2)/max(A(:,2))) ;
% disp(ratio.*ones(size(A(:,2)))) ;
% [xi,yi] = polyxpoly(A(:,2)/max(A(:,2)), A(:,5)./A(:,6),A(:,7)/max(A(:,7)), 0.5.*ones(size(A(:,7)/max(A(:,7)))));

% r0: 
%[xiArc,yiArc] = polyxpoly(A(:,2)/max(A(:,2)), A(:,5)./A(:,6), A(:,2)/max(A(:,2)), ratio.*ones(size(A(:,2))) );

% r1: 
[xiArc,yiArc] = polyxpoly(A(:,2)/max(A(:,2)), (A(:,5)-A(:,4))./(A(:,6)-A(:,4)), A(:,2)/max(A(:,2)), ratio.*ones(size(A(:,2))) );
% find cross section point [xi(1), yi(1)], [xi(2), yi(2)] for arctic seal. 

disp("xiArc")
disp(xiArc)
disp("yiArc")
disp(yiArc)

xArcPosi        = linspace(xiArc(1)*max(A(:,2)), xiArc(2)*max(A(:,2)), n); 
xArcPeri        = pchip(A(:,2), A(:,3), xArcPosi);
xArcAreaWoshunt = pchip(A(:,2), (A(:,5)-A(:,4)), xArcPosi);
xArcAreaWshunt  = pchip(A(:,2), (A(:,6)-A(:,4)), xArcPosi);

% r0: 
%[xiSub,yiSub] = polyxpoly(S(:,7)/max(S(:,7)), S(:,10)./S(:,11), S(:,7)/max(S(:,7)), ratio.*ones(size(S(:,7))));

% r1: 
[xiSub,yiSub] = polyxpoly(S(:,7)/max(S(:,7)), (S(:,10)-S(:,9))./(S(:,11)-S(:,9)), S(:,7)/max(S(:,7)), ratio.*ones(size(S(:,7))));
% find cross section point [xi(1), yi(1)], [xi(2), yi(2)] for subtropical seal. 

% disp("xiSub")
% disp(xiSub)
% disp("yiSub")
% disp(yiSub)


xSubPosi        = linspace(xiSub(1)*max(S(:,7)), xiSub(2)*max(S(:,7)), n); 
xSubPeri        = pchip(S(:,7), S(:,8),  xSubPosi);
xSubAreaWoshunt = pchip(S(:,7), (S(:,10)-S(:,9)), xSubPosi);
xSubAreaWshunt  = pchip(S(:,7), (S(:,11)-S(:,9)), xSubPosi);



%% Find the relation of dosal meatus area and specie  
%figure(2) 
%plot(A(:,2)/max(A(:,2)), A(:,6)-A(:,5), 'DisplayName','Arctic' )  ; hold on
%plot(S(:,7)/max(S(:,7)), S(:,11)-S(:,10), 'DisplayName','Subtropical' )  ;
%xlabel('z/L [-]','FontSize',14);
%ylabel('A_{D} [mm^2]','FontSize',14);
%legend('FontSize',14); hold off
% saveas(figure(2),'./fig/areaperi_subtract.pdf')
% saveas(figure(2),'./fig/areaperi_subtract.eps', 'epsc')


%figure(3)
%plot(A(:,2)/max(A(:,2)), A(:,6), 'DisplayName','Arctic' )  ; hold on
%plot(S(:,7)/max(S(:,7)), S(:,11), 'DisplayName','Subtropical' )  ; 
%%%plot(A(:,7)/A(end,7), A(:,11), 'DisplayName','Subtropical' )  ; 
%scatter(xArcPosi/max(A(:,2)), xArcAreaWshunt, 'DisplayName','Arctic' )  ; 
%scatter(xSubPosi/max(S(:,7)), xSubAreaWshunt, 'DisplayName','Subtropical' )  ; 
% disp("what is A(end,2) ?? ")
% disp(S(end,7))
% disp(S(:,7)./S(end,7))
% disp(S(:,11))
%xlabel('z/L [-]','FontSize',14);
%ylabel('A [mm^2]','FontSize',14);

%legend('Location', 'Best','FontSize',14); 
%saveas(figure(3),'./fig/area_subtract_new.pdf')
%saveas(figure(3),'./fig/area_subtract_new.eps', 'epsc')


%figure(4)
%plot(A(:,2)/max(A(:,2)), A(:,3), 'DisplayName','Arctic' )  ; hold on
%plot(S(:,7)/max(S(:,7)), S(:,8), 'DisplayName','Subtropical' )  ; 
%scatter(xArcPosi/max(A(:,2)), xArcPeri, 'DisplayName','Arctic' )  ; 
%scatter(xSubPosi/max(S(:,7)), xSubPeri, 'DisplayName','Subtropical' )  ; 
%xlabel('z/L [-]','FontSize',14);
%ylabel('\gamma [mm]','FontSize',14);

%legend('Location', 'northwest','FontSize',14);
%saveas(figure(4), [pwd '/fig/peri_subtract_new.pdf'])
%saveas(figure(4),'./fig/peri_subtract_new.eps', 'epsc')

result.xArc = xiArc;
result.xArcPosi = xArcPosi; 
result.xArcPeri = xArcPeri ;
result.xArcAreaWshunt = xArcAreaWshunt;
result.xArcAreaWoshunt = xArcAreaWoshunt;

result.xSub = xiSub; 
result.xSubPosi = xSubPosi; 
result.xSubPeri = xSubPeri ;
result.xSubAreaWshunt = xSubAreaWshunt;
result.xSubAreaWoshunt = xSubAreaWoshunt;

%result = [xArc; xArcPeri; xArcAreaWshunt; xSub; xSubPeri; xSubAreaWshunt];
end
