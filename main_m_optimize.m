clc, clear,close all;
%%
params.r = 0.5;
N = 20; % # of antennas
K = 8; %  # of users
exp_num = 10;
M = 20; 
R_set=1:0.5:8;
params.iter_max =10;
params.rho = 5;
params.E = 2;
params.verb = 1;
params.snr = 10^3; %30dB
iter_max = 10;
% [m,~, mse_sdr] = alterMin_SDR(Hd,Hr,G,iter_max,params);

m_IRS= nan(length(R_set),1);
 m_wo_IRS = nan(length(R_set),1);
m_rndV_IRS=nan(length(R_set),1);

for ii = 1:length(R_set)
    params.r=R_set(ii);
    tmp_m_IRS = 0;
   tmp_m_wo_IRS = 0;
   tmp_m_rndV_IRS=0;
 
   
    fprintf('Eset = %d\n',ii)
    parfor jj = 1:exp_num
%         fprintf('Nset = %d, exp_num = %d\n',ii,jj)
      Hd = normrnd(0,1/sqrt(2),N,K)+1i* normrnd(0,1/sqrt(2),N,K); %channel user to FC
      Hr = normrnd(0,1/sqrt(2),M,K)+1i* normrnd(0,1/sqrt(2),M,K); %channel user to IRS
      G  = normrnd(0,1/sqrt(2),N,M)+1i* normrnd(0,1/sqrt(2),N,M); %channe IRS to FC

%         
%         [~, mse_wo,~] = find_M_DC(Hd,params)
%         tmp_wo_IRS = tmp_wo_IRS+mse_wo;
     
        [m_wo_IRS,~,~] = find_minsum_m(Hd,params);
        tt=norm(m_wo_IRS);
       tmp_m_wo_IRS = tmp_m_wo_IRS+tt;
         
       [m_IRS,~,~] = alterminsum(Hd,Hr,G,iter_max,params)
            tmp1 = sum(~isnan(m_IRS));
%         tmp2 = sum_set2(~isnan(sum_set2));
%          tmp_m_rndV_IRS= tmp_m_rndV_IRS+tmp1(1);
        tmp_m_IRS =tmp_m_IRS+tmp1(end); 
        
    end
  m_IRS(ii) = tmp_m_IRS/exp_num;
  m_wo_IRS(ii)=tmp_m_wo_IRS/exp_num;
%    m_rndV_IRS (ii)= tmp_m_rndV_IRS/exp_num;
end
save m_optimize.mat 
%%
figure;
semilogy(R_set,m_IRS, '*-','LineWidth',2,'MarkerSize',12) 
hold on;
semilogy(R_set,m_wo_IRS, 'o-','LineWidth',2,'MarkerSize',12) 
hold on;
% semilogy(E_set,m_rndV_IRS, 'm-','LineWidth',2,'MarkerSize',12) 
% % hold on;
xlabel('RRR ','FontSize',14)
ylabel('b','FontSize',14)
xlim([1 8]);
% legend('sum power constraint','peak power constraint')
legend('Optimal reception beamforming b with IRS','Optimal reception beamforming b without IRS')
set(gca,'xtick',[1 2  3 4  5 6 8 ],'xticklabel',[1 2 3 4 5 6 8])
hold off

