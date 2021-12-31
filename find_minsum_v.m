function [y,feasibility] = find_minsum_v(m,Hd,Hr,G,params)
iter_max = params.iter_max;
verb = params.verb;
[~,K] = size(Hd);
[~,M] = size(G);
feasibility =1;
tmp = randn(M+1,M+1)+1i*randn(M+1,M+1);
V = tmp*tmp';
[u,~] = eigs(V); %返回对角矩阵s 和矩阵 u，前者包含主对角线上的特征值，后者的各列中包含对应的特征向量。
V_partial = u(:,1)*u(:,1)';%特征值uu'
obj0 = 0;
for  iter = 1:iter_max
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
          % 1- power*(real(trace(R*V)+c'*c)) <= 0;
      tf=tf+inv_pos(real(trace(R*V)+c'*c));
       end
    minimize (tf+real(trace((eye(M+1)-V_partial')* V)))%I = eye(n) 返回一个主对角线元素为 1 且其他位置元素为 0 的 n×n 单位矩阵。
%     minimize (M+1-real(trace(V_partial'* V)))
     subject to      

        real(diag(V)) == 1;
    cvx_end
     
    if strcmp(cvx_status,'Infeasible') == 1 
        feasibility = 0;
        y = nan;
        return;
    end
    err = abs(cvx_optval-obj0);
    obj0 = cvx_optval;

    [u,s] = eigs(V);%返回对角矩阵s 和矩阵 u，前者包含主对角线上的特征值，后者的各列中包含对应的特征向量。
    V_partial = u(:,1)*u(:,1)';
    res = abs(norm(V,'fro')-s(1,1));  %norm(X,'fro') 返回矩阵 X 的 Frobenius 范数。
%     res = trace(s)-s(1,1);
    if verb>=2
        fprintf(' iter:%d/%d, err:%.3e, res:%.3e\n', iter, iter_max, err, res);
    end
    if err<1e-9 && res<1e-8 
        break;
    end
%     if  res<1e-8 
%         break;
%     end
end
[u,s,v] = svd(V);
y_bar = u(:,1)*sqrt(s(1,1));
y = y_bar(1:M)/y_bar(end);
% feasibility = check_feasible(y,m,Hd,Hr,G);
% y = y./abs(y);

end