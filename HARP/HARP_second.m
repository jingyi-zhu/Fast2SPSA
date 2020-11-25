function [theta_hat_ks, loss_ks] = ...
    HARP_second(a,A,alpha,c,gamma,loss,loss_noisy,n,rep,theta_0)

p = length(theta_0);
%% initializaiton
theta_hat_ks_all = zeros(p, n, rep);
loss_ks_all = zeros(1, n, rep);

eig_threshold = 1e-6;

%% iteration
for rep_idx = 1:rep
    H_bar_k = eye(p);
    H_barbar_k = eye(p);
    theta_hat_k = theta_0;
    for k = 0:(n-1)
        a_k = a / (k+1+A)^alpha; w_k = 1e-1/ (k+2);
        c_k = c / (k+1)^gamma; c_tilde_k = c / (k+1)^gamma;
        
        H_sqrt = sqrtm(H_barbar_k);  temp1 = 2 * round(rand(p, 1)) - 1; delta_k = H_sqrt \ temp1; tdelta_k = H_sqrt * temp1;
        %         delta_k = mvnrnd(zeros(p, 1), inv(H_barbar_k))'; tdelta_k = delta_k;
        temp2 = 2 * round(rand(p, 1)) - 1; delta_tilde_k = H_sqrt \ temp2; tdelta_tilde_k = H_sqrt * temp2;
        %         delta_tilde_k = mvnrnd(zeros(p, 1), inv(H_barbar_k))'; tdelta_tilde_k = delta_tilde_k;
        
        %%
        theta_hat_k_plus = theta_hat_k + c_k * delta_k; y_hat_k_plus = loss_noisy(theta_hat_k_plus);
        theta_hat_k_minus = theta_hat_k - c_k * delta_k; y_hat_k_minus = loss_noisy(theta_hat_k_minus);
        g_hat_k = (y_hat_k_plus - y_hat_k_minus) / (2 * c_k) .* tdelta_k;
        
        theta_hat_k = theta_hat_k - a_k * ( H_barbar_k \  g_hat_k );
        
        %%
        theta_hat_k_plus_perturbed = theta_hat_k_plus + c_tilde_k * delta_tilde_k; y_hat_k_plus_perturbed = loss_noisy(theta_hat_k_plus_perturbed);
        theta_hat_k_minus_perturbed = theta_hat_k_minus + c_tilde_k * delta_tilde_k; y_hat_k_minus_perturbed = loss_noisy(theta_hat_k_minus_perturbed);
        delta_y = (y_hat_k_plus_perturbed - y_hat_k_plus) - (y_hat_k_minus_perturbed - y_hat_k_minus);
        H_hat_k = (tdelta_tilde_k * tdelta_k'  + tdelta_k * tdelta_tilde_k') * delta_y / (4 * c_k * c_tilde_k);
        
        %%
        H_bar_k = (1 - w_k) * H_bar_k + w_k * H_hat_k;
        % H_barbar_k = sqrtm(H_bar_k' * H_bar_k + eig_threshold * eye(p));
        % [H_bar_k_V, H_bar_k_D] = eig(H_bar_k); H_bar_k_D = diag(H_bar_k_D); ind_neg = find(H_bar_k_D<=0); H_bar_k_D(ind_neg)  = max(eig_threshold, abs(H_bar_k_D(ind_neg)));  H_barbar_k = H_bar_k_V * diag(H_bar_k_D) * H_bar_k_V'; H_barbar_k = 0.5 * (H_barbar_k + H_barbar_k');
        [H_bar_k_V, H_bar_k_D] = eig(H_bar_k); H_bar_k_D = diag(H_bar_k_D); H_bar_k_D  = max(eig_threshold, abs(H_bar_k_D));  H_barbar_k = H_bar_k_V * diag(H_bar_k_D) * H_bar_k_V'; H_barbar_k = 0.5 * (H_barbar_k + H_barbar_k');
        
        %% record result
        theta_hat_ks_all(:, k+1, rep_idx) = theta_hat_k;
        loss_ks_all(1, k+1, rep_idx) = loss(theta_hat_k);
    end
end

theta_hat_ks = mean(theta_hat_ks_all, 3); % dim: p-by-n
loss_ks = mean(loss_ks_all, 3); % dim: 1-by-n
end