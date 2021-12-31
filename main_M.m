clc, clear,close all;
%%
% rng(0);
N = 10; % # of antennas
K = 8; %  # of users
exp_num = 10;
M_set = 2:2:12; 
params.r = 1;
params.iter_max =10;
params.rho = 1;
params.E = 0.2;
params.verb = 1;
params.snr = 10^4; %30dB
iter_max = 10;
%  [m,~, mse_sdr] = alterMin_SDR(Hd,Hr,G,iter_max,params);

sum_DCpower = nan(length(M_set),1);
 sum_wo_power = nan(length(M_set),1);
%  sum_rndV_power=nan(length(M_set),1);
 sum_SDR_power = nan(length(M_set),1);
for ii = 1:length(M_set)
    M=M_set(ii);
    tmp_DC_IRS = 0;
   tmp_wo_IRS = 0;
%    tmp_rndV_IRS=0;
   tmp_SDR_IRS = 0;
  
    fprintf('Mset = %d\n',ii)
    parfor jj = 1:exp_num
%         fprintf('Nset = %d, exp_num = %d\n',ii,jj)
      Hd = normrnd(0,1/sqrt(2),N,K)+1i* normrnd(0,1/sqrt(2),N,K); %channel user to FC
      Hr = normrnd(0,1/sqrt(2),M,K)+1i* normrnd(0,1/sqrt(2),M,K); %channel user to IRS
      G  = normrnd(0,1/sqrt(2),N,M)+1i* normrnd(0,1/sqrt(2),N,M); %channe IRS to FC

%  [Hr,G,Hd]= channel_realization_IRS(K,M,N);
     
        [~,sum_wu_IRS,~] = find_minsum_m(Hd,params);
         tmp_wo_IRS = tmp_wo_IRS+sum_wu_IRS;
         
         
       [~,~,sum] = alterminsum(Hd,Hr,G,iter_max,params)
            tmp1 = sum(~isnan(sum));
       
%           tmp_rndV_IRS= tmp_rndV_IRS+tmp1(1);
        tmp_DC_IRS =tmp_DC_IRS+tmp1(end); 
      [~,~,sum] =  alterMin_SDR(Hd,Hr,G,iter_max,params)
            tmp3 = sum(~isnan(sum));
        tmp_SDR_IRS =tmp_SDR_IRS+tmp3(end); 
        
    end
  sum_DCpower(ii) =tmp_DC_IRS/exp_num;
  sum_wo_power (ii)=tmp_wo_IRS/exp_num;
  sum_SDR_power (ii)=tmp_SDR_IRS/exp_num;
%      sum_rndV_power (ii)= tmp_rndV_IRS/exp_num;
end
save mainmmmm222.mat 
%%
figure;

semilogy(M_set,sum_wo_power, 'o-','LineWidth',2,'MarkerSize',12) 
hold on;
semilogy(M_set,sum_DCpower, '*-','LineWidth',2,'MarkerSize',12) 
hold on;
semilogy(M_set,sum_SDR_power, 'v-','LineWidth',2,'MarkerSize',12) 
 hold on;
% semilogy(M_set,sum_rndV_power, 'm-','LineWidth',2,'MarkerSize',12) 
%  hold on;
xlabel('number of elements at IRS M','FontSize',14)
ylabel('sumpower','FontSize',14)


legend('sum power without IRS','sum power DC','sum power SDR')

grid on


