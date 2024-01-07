clc
close all
clear all

lambda = 5;

w=0.0052*pi; 
A=[cos(w) -sin(w);sin(w) cos(w)]; 
H=[1 0];

x0 = [3;0];
P0 = 0.01.*eye(2);

Q = 0.01*eye(2);
R = 1;

num_test = 1e3;

mse_e_store = cell(1,num_test);
x_store = cell(1,num_test);
y_store = cell(1,num_test);

mse_e = zeros(2,1000);
mse_e_usm = mse_e;
mse_e_usm_upt = mse_e;
mse_e_pf = mse_e;
mse_e_conv = mse_e;
mse_e_tkf = mse_e;



for test = 1 : num_test
    
    if mod(test,10) == 0
        test
    end

for t = 1 : 1000
    if t == 1
        x(:,t) = A * x0 + sqrt(Q)* randn(2,1);
    else
        x(:,t) = A * x(:,t-1) + sqrt(Q)* randn(2,1);
    end
    
    y(:,t) = H*x(:,t) + sqrt(R) * randn(1,1);
    z(:,t) = modulo_adc(y(:,t),lambda);
    z_sat(:,t) = sat_adc(y(:,t));
    
end

y_store{1,test} = [y;z;z_sat];

% filtering

x0_est(:,test) = x0 + sqrt(P0)*randn(2,1);

% unlimited sensing filter

[x_e_usm,p_t] = kf_usm(A,Q,H,R,z,x0_est(:,test),P0,lambda);

% particle filter

x_e_pf = pf_usm(A,Q,H,R,z,x0_est(:,test),P0,lambda);


% conventional Kalman filter

x_e_conv = kf_conv(A,Q,H,R,z,x0_est(:,test),P0);


% tobik Kalman filter

x_e_tkf = kf_tobit(A,Q,H,R,z_sat,x0_est(:,test),P0);


mse_e_usm = mse_e_usm + (x - x_e_usm).^2;
mse_e_pf = mse_e_pf + (x - x_e_pf).^2;
mse_e_conv = mse_e_conv + (x - x_e_conv).^2;
mse_e_tkf = mse_e_tkf +(x - x_e_tkf).^2;

x_store{1,test} = [x;x_e_usm;x_e_pf;x_e_conv;x_e_tkf];
mse_e_store{1,test} = [(x - x_e_usm).^2;(x - x_e_pf).^2;...
                          (x - x_e_conv).^2;(x - x_e_tkf).^2;];
end


figure
plot((mse_e_usm(1,:)./num_test))
hold on
plot((mse_e_pf(1,:)./num_test))
hold on
plot((mse_e_conv(1,:)./num_test))
hold on
plot((mse_e_tkf(1,:)./num_test))



figure
plot((mse_e_usm(2,:)./num_test))
hold on
plot((mse_e_pf(2,:)./num_test))
hold on
plot((mse_e_conv(2,:)./num_test))
hold on
plot((mse_e_tkf(2,:)./num_test))


save simulation_data.mat


    
    
    


