function [m,sum,feasibility] = find_minsum_m(H,params)
E=params.E;
r=params.r;
maxiter = params.iter_max;
verb = params.verb;
feasibility=1;
[N,K] = size(H);
tmp = randn(N,N)+1i*randn(N,N);%
M = tmp*tmp';
[u,~,~] = svd(M);
M_partial = u(:,1)*u(:,1)';%
obj0 = 0;

for  iter = 1:maxiter
    %% Solve convex subproblem
    cvx_begin quiet
    cvx_solver sdpt3
    variable M(N,N) hermitian semidefinite
      expression tt          
      tt=0;
          for k=1:K
              tt=tt+inv_pos(real(H(:,k)'*M*H(:,k))) ;
         % 1-(real(H(:,k)'*M*H(:,k)) )<= 0 ;
          end
%     minimize((1+rho)*trace(M)- real(trace(rho*M_partial'* M)))
    minimize (tt+real(trace((eye(N)- M_partial')* M)))
     subject to  

          real(trace(M))<= E ;
real(trace(M))<= 1/(2^r-1) ;  
    cvx_end
    
    if strcmp(cvx_status,'Infeasible')
        feasibility = 0;
        m = nan;
        return;
    end
    
    err = abs(cvx_optval-obj0);
    obj0 = cvx_optval;
    
    %% Subgradient
    [u,s] = eigs(M,1);
    M_partial = u(:,1)*u(:,1)';
    res = abs(norm(M,'fro')-s(1,1));
%     res = trace(s)-s(1,1);
    if verb>=2
        fprintf(' iter:%d/%d, err:%.3e, res:%.3e\n', iter, maxiter, err, res);
    end
    if err<1e-13
        break;
    end
    if  res<1e-8 && err<1e-8
        break;
    end
end
[u,s,v] = svd(M);
m = u(:,1)*sqrt(s(1,1));

 
%% compute sum power
sum = 0;
for iter=1:K
    sum =sum+ 1/norm(m'*H(:,iter))^2;
  
end

end