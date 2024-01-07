function [x_e] = kf_conv(A,Q,H,R,z,x0_est,P0)

for t = 1 : 1000
    
    if t == 1
        x_pred(:,t) = A*x0_est;
        P_pred = A * P0 * A' + Q;
    else
        x_pred(:,t) = A*x_e(:,t-1);
        P_pred = A * P_e * A' + Q;
    end
    
    K=P_pred*H'/(H*P_pred*H'+R); 
    
    P_e=(eye(2)-K*H)*P_pred;
    
    y1_pred = H*x_pred(:,t);
    x_e(:,t) = x_pred(:,t) + K*(z(:,t)  - y1_pred);
    
end

end