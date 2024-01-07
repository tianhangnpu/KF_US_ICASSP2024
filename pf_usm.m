% Particle filter

function xf = pf_usm(A,Q,H,R,z,x0_est,P0,lambda)

N = 1000;
M = 1000;  %number of samples
xf = zeros(2,1000);
xPWeighted  = zeros(N,M);      % Make some room


x = repmat(x0_est,1,M) + sqrt(P0)*randn(2,M);

for t=1:N
    
    for i = 1 : M
        e_ff(:,i) = modulo_adc(H*x(:,i),lambda);
    end
    
    
    e = repmat(z(t),1,M) - e_ff;
    
    
    w = exp(-(1/2)*(e.*(R\e)));  % 2. Evaluate the likelihood  % m.R\e = inv(m.R)*e %
    w = w/sum(w);
    
    xf(:,t) = sum(repmat(w,2,1).*x,2);             % Compute state estimate
    index = sysresample(w);        % 3. Resample
    x     = x(:,index);
    
    x = A * x + sqrt(Q)*randn(2,M); % 4. Predict the particles
end;



end


