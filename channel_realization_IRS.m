function [Hr,G,Hd]= channel_realization_IRS(K,M,N)
%K: number of user
%M: number of elements at IRS
%N: number of antenna at AP

location_AP = [0,0,0];
location_IRS = [20,20,20];

%path loss parameter
T0=100000;
a_AI = 2.2;
a_AU = 3.5;
a_IU = 2.8;

dis_AI = norm(location_AP-location_IRS);
Ploss_AI = path_loss(T0,dis_AI,1,a_AI);

Ploss_AU = nan(K,1);
Ploss_IU = nan(K,1);
for k = 1:K
    x = 30 + 60*rand;
    y = -30+ 60*rand;
    location_user = [x,y,0];
    
    dis_AU = norm(location_AP-location_user);
    Ploss_AU(k) = path_loss(T0,dis_AU,1,a_AU);
    
    dis_IU = norm(location_IRS-location_user);
    Ploss_IU(k) = path_loss(T0,dis_IU,1,a_IU);
end

%% small-scale fading
G =  normrnd(0,1/sqrt(2),N,M)+1i* normrnd(0,1/sqrt(2),N,M);
G = sqrt(Ploss_AI)*G; %channe IRS to AP

Hd = normrnd(0,1/sqrt(2),N,K)+1i* normrnd(0,1/sqrt(2),N,K);
Hr = normrnd(0,1/sqrt(2),M,K)+1i* normrnd(0,1/sqrt(2),M,K); 
for k = 1:K
Hd(:,k) = sqrt(Ploss_AU(k)) * Hd(:,k);  %channel user to AP
Hr(:,k) = sqrt(Ploss_IU(k)) * Hr(:,k); %channel user to IRS
end

end

function Ploss = path_loss(T0,d,d0,a)
Ploss = T0*(d/d0)^(-a);
end