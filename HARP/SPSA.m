function [theta_hat_ks, loss_ks] = ...
    SPSA(a,A,alpha,c,gamma,loss,loss_noisy,n,rep,theta_0)

p = length(theta_0);

%% initializaiton
theta_hat_ks_all = zeros(p, n, rep);
loss_ks_all = zeros(1, n, rep);

%% iteration
for rep_idx = 1:rep
    for k = 0:(n-1)
        if k == 0
            theta_hat_k = theta_0;
        end
        a_k = a / (k+1+A)^alpha; c_k = c / (k+1)^gamma;
        
        g_hat_k = zeros(p, 1);
        avg_num = 2;
        for avg = 1 : avg_num
            delta_k = 2 * round(rand(p, 1)) - 1;
            
            theta_hat_k_plus = theta_hat_k + c_k * delta_k;
            theta_hat_k_minus = theta_hat_k - c_k * delta_k;
            
            y_hat_k_plus = loss_noisy(theta_hat_k_plus);
            y_hat_k_minus = loss_noisy(theta_hat_k_minus);
            g_hat_k = g_hat_k + (y_hat_k_plus - y_hat_k_minus) ./ (2 * c_k * delta_k) / avg_num;
        end
        
        theta_hat_k = theta_hat_k - a_k * g_hat_k;  
        
        theta_hat_ks_all(:, k+1, rep_idx) = theta_hat_k;
        loss_ks_all(1, k+1, rep_idx) = loss(theta_hat_k);
    end
end

theta_hat_ks = mean(theta_hat_ks_all, 3); % dim: p-by-n
loss_ks = mean(loss_ks_all, 3); % dim: 1-by-n
end