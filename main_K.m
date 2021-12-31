clc, clear,close all;
%%
% rng(0);
N = 10; % # of antennas
K_set = 2:1:16; %  # of users
exp_num = 10;
M = 20; 

params.r = 0.5;
params.iter_max =10;
params.rho = 5;
params.E = 2;
params.verb = 1;
params.snr = 10^3; %30dB
iter_max = 10;
% [m,~, mse_sdr] = alterMin_SDR(Hd,Hr,G,iter_max,params);

sum_power = nan(length(K_set),1);
 sum_wo_power = nan(length(K_set),1);
sum_rndV_power=nan(length(K_set),1);

for ii = 1:length(K_set)
    K=K_set(ii);
    tmp_power_IRS = 0;
   tmp_wo_IRS = 0;
   tmp_rndV_IRS=0;
  
    fprintf('Kset = %d\n',ii)
    parfor jj = 1:exp_num
%         fprintf('Nset = %d, exp_num = %d\n',ii,jj)
      Hd = normrnd(0,1/sqrt(2),N,K)+1i* normrnd(0,1/sqrt(2),N,K); %channel user to FC
      Hr = normrnd(0,1/sqrt(2),M,K)+1i* normrnd(0,1/sqrt(2),M,K); %channel user to IRS
      G  = normrnd(0,1/sqrt(2),N,M)+1i* normrnd(0,1/sqrt(2),N,M); %channe IRS to FC

%         
%         [~, mse_wo,~] = find_M_DC(Hd,params)
%         tmp_wo_IRS = tmp_wo_IRS+mse_wo;
     
        [~,sum_wu_IRS,~] = find_minsum_m(Hd,params);
         tmp_wo_IRS = tmp_wo_IRS+sum_wu_IRS;
         
       [~,~,sum] = alterminsum(Hd,Hr,G,iter_max,params)
            tmp1 = sum(~isnan(sum));
%         tmp2 = sum_set2(~isnan(sum_set2));
%          tmp_rndV_IRS= tmp_rndV_IRS+tmp1(1);
        tmp_power_IRS =tmp_power_IRS+tmp1(end); 
        
    end
  sum_power(ii) =tmp_power_IRS/exp_num;
  sum_wo_power (ii)=tmp_wo_IRS/exp_num;
%    sum_rndV_power (ii)= tmp_rndV_IRS/exp_num;
end
save sumpowerkkkkkkk22222.mat 
%%
figure;
semilogy(K_set,sum_power, '*-','LineWidth',2,'MarkerSize',12) 
hold on;
semilogy(K_set,sum_wo_power, 'o-','LineWidth',2,'MarkerSize',12) 
hold on;
% semilogy(K_set,sum_rndV_power, 'm-','LineWidth',2,'MarkerSize',12) 
% hold on;
xlim([2 16]);
xlabel('Number of wireless sensors K ','FontSize',14)
ylabel('sumpower','FontSize',14)

% legend('sum power constraint','peak power constraint')
legend('sum power with IRS','sum power without IRS','sum power with random IRS')
legend('sum power with IRS','sum power without IRS')
set(gca,'xtick',[2 4  6 8  10 12 14 16 ],'xticklabel',[2 4  6 8  10 12 14 16 ])
grid on


