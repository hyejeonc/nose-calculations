function sol = DATA_Residual(fileName, n, L, Nt, nFirst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% === input === 
% fileName: str, file directory 
% n       : int, discretization number for plotting
% L       : float, length of a nose 
% Nt      : int, discretization number for time
% 
% === output ===
% sol.(variable below)
%
% totSigma: array, 1*Ns 
% inSigma : array, 1*Ns
% outSigma: array, 1*Ns
% Ta025   : array, 1*Ns
% Ta050   : array, 1*Ns
% Ta075   : array, 1*Ns
% Ta100   : array, 1*Ns
% 
%
% plotTotSigma: array, 1*n 
% plotInSigma : array, 1*n
% plotOutSigma: array, 1*n
% plotTa025   : array, 1*n
% plotTa050   : array, 1*n
% plotTa075   : array, 1*n
% plotTa100   : array, 1*n
% 
%
% timeSigma     : array, 1*NtMax  
% timeTa075     : array, 1*NtMax
% timeResSigma  : array, 1*NtMax  
% timeResTa075  : array, 1*NtMax
% timeResTaNorm : array, 1*NtMax
% Ncyc: int  
%
% timeResDynaSigmaExp
% timeResDynaTaExp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    nose = load(fileName,'Nose'); 
    Ns = length(nose.Nose.x);
    NtMax = length(nose.Nose.t) / (2*Nt);
catch
    nose = load(fileName,'Sol'); 
    nose.Nose = nose.Sol;
    Ns = length(nose.Sol.x);
    NtMax = length(nose.Sol.t) / (2*Nt);
end
    
    

%Nt = length(); 


sol.Ncyc = NtMax;

for j = 1:Ns
    sol.totSigma(j) = trapz(nose.Nose.t(end-Nt*2+1:end),    nose.Nose.sigma(end-Nt*2+1:end,j)); 
    sol.inSigma(j)  = trapz(nose.Nose.t(end-Nt*2+1:end-Nt), nose.Nose.sigma(end-Nt*2+1:end-Nt,j)); 
    sol.outSigma(j) = trapz(nose.Nose.t(end-Nt+1:end)     , nose.Nose.sigma(end-Nt+1:end,j)); 

    
    sol.totSigmaHalf(j) = trapz(nose.Nose.t(end-Nt*2+1-(2*Nt)*NtMax/2:end),    nose.Nose.sigma(end-Nt*2+1-(2*Nt)*NtMax/2:end,j)); 
    sol.totSigmaQuar(j) = trapz(nose.Nose.t(end-Nt*2+1-(2*Nt)*NtMax/4:end),    nose.Nose.sigma(end-Nt*2+1-(2*Nt)*NtMax/4:end,j)); 

    
    sol.Ta025(j) = nose.Nose.Ta(end-Nt*(1.5),j);
    sol.Ta050(j) = nose.Nose.Ta(end-Nt*(1)  ,j);
    sol.Ta075(j) = nose.Nose.Ta(end-Nt*(0.5),j);
    sol.Ta100(j) = nose.Nose.Ta(end         ,j);
    
    sol.Ta025Half(j) = nose.Nose.Ta(end-Nt*(1.5)-(2*Nt)*NtMax/2,j);
    sol.Ta075Half(j) = nose.Nose.Ta(end-Nt*(0.5)-(2*Nt)*NtMax/2,j);
    
    sol.Ta025Quar(j) = nose.Nose.Ta(end-Nt*(1.5)-(2*Nt)*NtMax/4,j);
    sol.Ta075Quar(j) = nose.Nose.Ta(end-Nt*(0.5)-(2*Nt)*NtMax/4,j);  
end

%disp(sol.totSigma) ;% 12
disp("sum of tot Sigma ") ;% 200

disp(sum(sol.totSigma)) ; % 12
disp(sum(sol.outSigma) + sum(sol.inSigma));


for i = 1:NtMax
    
    % calculate total entropy for every cycle. 
    dsum = 0;
    for j = 1:Ns
%         disp("i")
%         disp(i)
%         disp(1+(i-1)*(2*Nt))
%         disp(i*(2*Nt))
    
        dsum = dsum + trapz(nose.Nose.t(1+(i-1)*(2*Nt):i*(2*Nt)), nose.Nose.sigma(1+(i-1)*(2*Nt):i*(2*Nt),j)); 

    end
    sol.timeSigma(i) = dsum;
    

    % calculate residual of total entropy for every cycle, non-scaled.
    if i > 1 
        sol.timeResSigma(i-1) = abs(sol.timeSigma(i) - sol.timeSigma(i-1)) ;
    end

    % calculate expelled Temperature of air at 0.75t for every cycle. x = 1
    %sol.timeTa025(i) = nose.Nose.Ta(0.25*(2*Nt)+(i-1)*(2*Nt),1) ;
    %sol.timeTa050(i) = nose.Nose.Ta(0.50*(2*Nt)+(i-1)*(2*Nt),1) ;
    sol.timeTa075(i) = nose.Nose.Ta(0.75*(2*Nt)+(i-1)*(2*Nt),1) ;
    %sol.timeTa100(i) = nose.Nose.Ta(     (2*Nt)+(i-1)*(2*Nt),1) ;
    
    % calculate residual of expelled Temperature of air at 0.75 for every cycle, non-scaled.
    if i > 1  
        sol.timeResTa075(i-1) = abs(sol.timeTa075(i) - sol.timeTa075(i-1)) ; 
    end
    
    
    
    % calculate residual of Temperature of air norm for every cycle, non-scaled.
    % sol.timeTaNorm(i) = norm(nose.Nose.Ta(1+(i-1)*(2*Nt):i*(2*Nt),:) - nose.Nose.Ta(1+(i)*(2*Nt):(i+1)*(2*Nt),:)) ;
    if i > 1  
        sol.timeResTaNorm(i-1) = norm(nose.Nose.Ta(1+(i-1)*(2*Nt):i*(2*Nt),:) - nose.Nose.Ta(1+(i-2)*(2*Nt):(i-1)*(2*Nt),:));

        %sol.timeResTaNorm(i-1) = sol.timeTaNorm(i-1) ; 
    end
   
end
disp("what is  tot Sigma at end, end-1 ?  ") % 200
disp(sol.timeSigma(end)) % 12
disp(sol.timeSigma(end-1))

for i = 1:nFirst
    for j = 1:2*Nt
        sol.timeResDynaSigmaExp(j+(i-1)*2*Nt) = sum(nose.Nose.sigma(j+(i-1)*Nt,:)); 
        sol.timeResDynaTaExp(j+(i-1)*2*Nt)    = nose.Nose.Ta(j+(i-1)*Nt,1);
    end
end


z = linspace(0,L,n);

sol.plotTotSigma = pchip(nose.Nose.x, sol.totSigma, z);
sol.plotTotSigmaQuar = pchip(nose.Nose.x, sol.totSigmaQuar, z);
sol.plotTotSigmaHalf = pchip(nose.Nose.x, sol.totSigmaHalf, z);

sol.plotInSigma  = pchip(nose.Nose.x, sol.inSigma, z);
sol.plotOutSigma = pchip(nose.Nose.x, sol.outSigma, z);
sol.plotTa025    = pchip(nose.Nose.x, sol.Ta025, z);
sol.plotTa050    = pchip(nose.Nose.x, sol.Ta050, z);
sol.plotTa075    = pchip(nose.Nose.x, sol.Ta075, z);
sol.plotTa100    = pchip(nose.Nose.x, sol.Ta100, z);

sol.plotTa025Half    = pchip(nose.Nose.x, sol.Ta025Half, z);
sol.plotTa075Half    = pchip(nose.Nose.x, sol.Ta075Half, z);
sol.plotTa025Quar    = pchip(nose.Nose.x, sol.Ta025Quar, z);
sol.plotTa075Quar    = pchip(nose.Nose.x, sol.Ta075Quar, z);

sol.timeResDynaSigmaExp
sol.timeResDynaTaExp
