function [y,feasibility ]= find_V_SDR(m,Hd,Hr,G,params)
[N,K] = size(Hd);
[~,M] = size(G);
E=params.E;
feasibility =1;
tmp = randn(M+1,M+1)+1i*randn(M+1,M+1);
V = tmp*tmp';

cvx_begin quiet
%                 cvx_solver sedumi
variable V(M+1,M+1) hermitian semidefinite
 expression tf
%     minimize (trace(V)- real(trace(V_partial'* V)))   
         tf=0;
       for k=1:K
            a_H = m'*G*diag(Hr(:,k));
            a = a_H';
            c = m'*Hd(:,k);
            R = [a*a',a*c;c'*a',0];
     
      tf=tf+inv_pos(real(trace(R*V)+c'*c));
       end
minimize (tf)
subject to

real(diag(V)) == 1;
cvx_end

if strcmp(cvx_status,'Infeasible') == 1 
    feasibility = 0;
    y = [];
    return; 
end

if rank(V,1e-6) == 1
    [u,s,v] = svd(V);
    y_bar = u(:,1)*sqrt(s(1,1));
    y = y_bar(1:M)/y_bar(end);
%     y = y./abs(y);
%     feasibility = check_feasible(y,m,Hd,Hr,G);
else
    for ii = 1:1000
        zi = chol(V,'lower');
        xi = (randn(M+1,1)+1i*randn(M+1,1))/sqrt(2);
        xi = zi*xi;        
       
%  y_bar = xi;
            y_bar = xi/E;
        end
        y = y_bar(1:M)/y_bar(end);
        y = y./abs(y);


    end       
end
