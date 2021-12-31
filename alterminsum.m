function [m_set,Theta,sum_set] = alterminsum(Hd,Hr,G,iter_max,params,theta)
[N,K] = size(Hd);
[~,M] = size(G);
E= params.E;
% sumpower2=0;
% snr = params.snr;
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
    
    if params.verb>=2
       fprintf('====iter: %d===find m =======\n',ii) 
    end

%     [m, ~] = find_minsum_m(He,params);
       [m,sum, ~] = find_minsum_m(He,params);
              
         sum_set(ii) = sum;  
       m_set(ii) = norm(m); 
    if ii>1 && abs(sum_set(ii)-sum_set(ii-1))<1e-3
        break;
    end
    
    if params.verb>=2
        fprintf('===iter: %d====find V =======\n',ii)
    end

 [Theta,isfeasible1] = find_minsum_v(m,Hd,Hr,G,params);
    if isfeasible1 ==0    
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