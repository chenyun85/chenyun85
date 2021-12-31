function [m_set,Theta,sum_set] = alterMin_SDR(Hd,Hr,G,iter_max,params,theta)
[N,K] = size(Hd);
[~,M] = size(G);

if nargin<6
    theta = randn(M,1)+1i*rand(M,1);
    theta = theta./abs(theta);
end
Theta = diag(theta);
He = nan(size(Hd));
for k=1:K
    He(:,k) = G*Theta*Hr(:,k)+Hd(:,k);
end
sum_set = nan(iter_max,1);
m_set= nan(iter_max,1);
for ii = 1:iter_max
    

    [m,sum,~] = find_M_SDR(He,params);
             sum_set(ii) = sum;  
       m_set(ii) = norm(m);  
    
    if ii>1 && abs(sum_set(ii)-sum_set(ii-1))<1e-3
        break;
    end
%     
%     if check_feasible(theta,m,Hd,Hr,G) == 0
%         error('The solution m is infeasible')
%     end
    
    [theta,isfeasible] = find_V_SDR(m,Hd,Hr,G,params);
    
    if isfeasible == 0
        break;
    else
        Theta = diag(theta);
        He = nan(size(Hd));
        for k=1:K
            He(:,k) = G*Theta*Hr(:,k)+Hd(:,k);
        end
    end
    
end
end