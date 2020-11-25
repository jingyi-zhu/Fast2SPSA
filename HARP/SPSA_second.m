function [theta_hat_ks, loss_ks] = ...
    SPSA_second(a,A,alpha,c,gamma,loss,loss_noisy,n,rep,theta_0) 

p = length(theta_0);

%% initializaiton
theta_hat_ks_all = zeros(p, n, rep);
loss_ks_all = zeros(1, n, rep);

%% iteration
for rep_idx = 1:rep
    Hbar_k = eye(p) * 1;
    %Hbar_k = diag([100; 1]);
    for k = 0:(n-1)
        if k == 0
            theta_hat_k = theta_0;
        end
        a_k = a / (k+1+A)^alpha;
        c_k = c / (k+1)^gamma;
        c_tilde_k = c / (k+1)^gamma;
        w_k = 1 / (k+2); % weight for new observation
        
        %%
        delta_k = 2 * round(rand(p, 1)) - 1;
        delta_tilde_k = 2 * round(rand(p, 1)) - 1;
        
        theta_hat_k_plus = theta_hat_k + c_k * delta_k;
        y_hat_k_plus = loss_noisy(theta_hat_k_plus);
        theta_hat_k_minus = theta_hat_k - c_k * delta_k;
        y_hat_k_minus = loss_noisy(theta_hat_k_minus);
        
        g_hat_k = (y_hat_k_plus - y_hat_k_minus) ./ (2 * c_k * delta_k);
        
        if k <= 0.2 * n - 1
            theta_hat_k = theta_hat_k - a_k * (g_hat_k);
        else
            theta_hat_k_plus_perturbed = theta_hat_k_plus + c_tilde_k * delta_tilde_k;
            y_hat_k_plus_perturbed = loss_noisy(theta_hat_k_plus_perturbed);
            
            theta_hat_k_minus_perturbed = theta_hat_k_minus + c_tilde_k * delta_tilde_k;
            y_hat_k_minus_perturbed = loss_noisy(theta_hat_k_minus_perturbed);
            G_hat_k_plus = (y_hat_k_plus_perturbed - y_hat_k_plus) ./ (c_tilde_k * delta_tilde_k);
            G_hat_k_minus = (y_hat_k_minus_perturbed - y_hat_k_minus) ./ (c_tilde_k * delta_tilde_k);
            delta_G_k = G_hat_k_plus - G_hat_k_minus;
            
            H_hat_k = delta_G_k' ./ (2 * c_k * delta_k);
            H_hat_k = (H_hat_k + H_hat_k') / 2;
            
            Hbar_k = (1 - w_k) * Hbar_k + w_k * H_hat_k;
            Hbarbar_k = sqrtm(Hbar_k' * Hbar_k + 1e-5 * eye(p));
            
            theta_hat_k = theta_hat_k - a_k * (Hbarbar_k \ g_hat_k);
            
%             %% blocking
%             temp = theta_hat_k - a_k * (Hbarbar_k \ g_hat_k);
%             if norm(temp - theta_hat_k) <= 2
%                 theta_hat_k = temp;
%             end
%             
%             %% projection
%             theta_hat_k = median([ -2.048 * ones(p, 1), theta_hat_k, 2.047 * ones(p, 1)], 2);
        end
        
        % record result
        theta_hat_ks_all(:, k+1, rep_idx) = theta_hat_k;
        loss_ks_all(1, k+1, rep_idx) = loss(theta_hat_k);
    end
end

theta_hat_ks = mean(theta_hat_ks_all, 3); % dim: p-by-n
loss_ks = mean(loss_ks_all, 3); % dim: 1-by-n
end