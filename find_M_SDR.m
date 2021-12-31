function [m, sum ,feasibility] = find_M_SDR(H,params)
[N,K] = size(H);
E=params.E;
r=params.r;
cvx_begin quiet
cvx_solver sdpt3
variable M(N,N) hermitian semidefinite
 expression tt          
      tt=0;
          for k=1:K
              tt=tt+inv_pos(real(H(:,k)'*M*H(:,k))) ;
          end
minimize (tt)
subject to
for k=1:K
            real(trace(M))<= E ;
real(trace(M))<= 1/(2^r-1) ; 
end
cvx_end
[u,s] = eigs(M);
m = u(:,1)*sqrt(s(1,1));
% m = u(:,1);
%% 
if rank(M,1e-6)>1
    zi = chol(M,'lower');
    xi = (randn(N,1)+1i*randn(N,1))/sqrt(2);
    xi = zi*xi;
    min_eig = inf;
    for k=1:K
        tmp = abs(xi'*xi)/E;
%          tmp = abs(xi);
        if tmp<min_eig
            min_eig = tmp;
        end
    end
    m = xi/min_eig;   
end

%% check feasibility
feasibility = 1;
% if check_feasible([],m,H) == 0
%     warning('The solution m is infeasible')
%     feasibility = 0;
%     mse = nan;
%     return;
% end
%  
%% compute MSE
sum = 0;
for iter=1:K
    sum =sum+ 1/norm(m'*H(:,iter))^2;
  
end


end