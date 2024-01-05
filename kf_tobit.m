function x_e = kf_tobit(A,Q,H,R,z,x0_est,P0)

sigma = sqrt(R);
tau = 0;

for t = 1 : 1000
    if t == 1
        x_pred(:,t) = A*x0_est;
        P_pred = A * P0 * A' + Q;
    else
        x_pred(:,t) = A*x_e(:,t-1);
        P_pred = A * P_e * A' + Q;
    end
    
    Ep=cdf('Normal',H*x_pred(:,t),tau,sigma);
    
    Rxy=P_pred*H'*Ep;
    
    pd=pdf('Normal',tau,H*x_pred(:,t),sigma);
    
    cd=cdf('Normal',tau,H*x_pred(:,t),sigma);
    
    temp = 1 - cd;
    
    if temp == 0
        temp = 1e-5;
    end
    
    m=pd/(temp);
    % temp = m * (m - (tau -(H*x_pred(:,t))/sigma ));
    temp = m * (m - ((tau -H*x_pred(:,t))/sigma ));
    
    Ryy=Ep*H*P_pred*H'*Ep+sigma*sigma*(1-(m*(m-((tau-H*x_pred(:,t))/sigma))));
    % Ryy=Ep*H*P_pred*H'*Ep+temp;
    Ey=Ep*(H*x_pred(:,t)+sigma*(m))+pdf('Normal',tau,H*x_pred(:,t),sigma)*tau;
    
    
    
    x_e(:,t)=x_pred(:,t)+Rxy/Ryy*(z(:,t)-Ey);
    P_e=(eye(2)-Ep*Rxy/Ryy*H)*P_pred;
end



end