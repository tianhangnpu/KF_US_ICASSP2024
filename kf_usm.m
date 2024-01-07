function [x_e,p] = kf_usm(A,Q,H,R,z,x0_est,P0,lambda)

offset = [-2,-1,0,1,2];

for t = 1 : 1000
    
    if t == 1
        x_pred(:,t) = A*x0_est;
        P_pred = A * P0 * A' + Q;
    else
        x_pred(:,t) = A*x_e(:,t-1);
        P_pred = A * P_e * A' + Q;
    end
    
    K=P_pred*H'/(H*P_pred*H'+R);
    
    S = (H*P_pred*H'+R);
    
    P_e=(eye(2)-K*H)*P_pred;
    
    y1_pred = H*x_pred(:,t);
    
    e_est = z(:,t) - y1_pred;
    
    e_est = e_est ./ lambda;
    
    e_est = ceil(e_est);
    
    e_possible = e_est + offset;

    
    x_possible = cell(1,5);
    
    for i = 1 : 5
        x_possible{1,i} = x_pred(:,t) + K*(z(:,t) - e_possible(i).*5 - y1_pred);
        p_temp(i,:) = exp(- (z(:,t) - e_possible(i).*5 - y1_pred).^2./S./2)./sqrt(2.*pi .* S);
    end
    
   
    
    p(:,t) = p_temp ./ sum(p_temp);
    x_temp = zeros(2,1);
    
    P_e_upt = zeros(size(P_e));
    
    p_t = p(:,t);
    
   
    
    for i = 1 : 5
        x_temp = x_temp +  p_t(i) .* x_possible{1,i};
        P_e_upt = P_e_upt + p_t(i) .* P_e ;
    end
    x_e(:,t) = x_temp;
    P_e = P_e_upt;
    
    
    
    
end

end