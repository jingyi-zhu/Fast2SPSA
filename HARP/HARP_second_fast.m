function [theta_hat_ks, loss_ks] = ...
    HARP_second_fast(a,A,alpha,c,gamma,loss,loss_noisy,n,rep,theta_0)
p = length(theta_0);
%% initializaiton
theta_hat_ks_all = zeros(p, n, rep);
loss_ks_all = zeros(1, n, rep);

eig_LB = 1e-6;

%% iteration
for rep_idx = 1:rep
    % P * H * P' = L * B * L'
    % P: permutation matrix
    % H: estimated Hession
    % L: lower triangle matrix
    % B: block diagonal matrix of size 2
    P_k_vec = 1:p; % store Q in vector form, index position of 1
    H_bar_k_factor = eye(p); % store M and B
    change_k = ones(1,p); % block structure indicator for B
    
    L_k = eye(p);
    % B_k_sqrt = V * D * V'
    B_k_sqrt_V = eye(p);
    B_k_sqrt_D_vec = ones(p,1);
    
    theta_hat_k = theta_0;
    
    for k = 0:(n-1)
        a_k = a / (k+1+A)^alpha; w_k = 1e-1/ (k+2); d_k = 1 - w_k;
        c_k = c / (k+1)^gamma; c_tilde_k = c / (k+1)^gamma;
        
        temp1 = 2 * round(rand(p, 1)) - 1;
        temp2 = 2 * round(rand(p, 1)) - 1;
        
        % compute H_sqrt \ temp
        delta_k = get_delta(P_k_vec, L_k, B_k_sqrt_V, B_k_sqrt_D_vec, temp1, p);
        delta_tilde_k = get_delta(P_k_vec, L_k, B_k_sqrt_V, B_k_sqrt_D_vec, temp2, p);
        % compute H_sqrt * temp
        tdelta_k = get_tdelta(P_k_vec, L_k, B_k_sqrt_V, B_k_sqrt_D_vec, temp1, p);
        tdelta_tilde_k = get_tdelta(P_k_vec, L_k, B_k_sqrt_V, B_k_sqrt_D_vec, temp2, p);
        
        %%
        theta_hat_k_plus = theta_hat_k + c_k * delta_k; y_hat_k_plus = loss_noisy(theta_hat_k_plus);
        theta_hat_k_minus = theta_hat_k - c_k * delta_k; y_hat_k_minus = loss_noisy(theta_hat_k_minus);
        g_hat_k = (y_hat_k_plus - y_hat_k_minus) / (2 * c_k) .* tdelta_k;
        
        theta_hat_k = theta_hat_k - a_k * ( H_barbar_k \ g_hat_k ) ;
                
        %%
        theta_hat_k_plus_perturbed = theta_hat_k_plus + c_tilde_k * delta_tilde_k; y_hat_k_plus_perturbed = loss_noisy(theta_hat_k_plus_perturbed);
        theta_hat_k_minus_perturbed = theta_hat_k_minus + c_tilde_k * delta_tilde_k; y_hat_k_minus_perturbed = loss_noisy(theta_hat_k_minus_perturbed);
        delta_y_k = (y_hat_k_plus_perturbed - y_hat_k_plus) - (y_hat_k_minus_perturbed - y_hat_k_minus);
        
        %%
        b_k = w_k * delta_y_k / (4 * c_k * c_tilde_k);
        u_k = tdelta_tilde_k; v_k = tdelta_k;
        u_k_norm = norm(u_k); v_k_norm = norm(v_k);
        u_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k + u_k_norm/v_k_norm*v_k);
        v_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k - u_k_norm/v_k_norm*v_k);
        
        % update H_bar_factor
        H_bar_k_factor(1:(p+1):(p^2)) = d_k * H_bar_k_factor(1:(p+1):(p^2));
        block_k_idx = find(change_k == 2);
        H_bar_k_factor(1+(block_k_idx-1)*(p+1)+1) = d_k * H_bar_k_factor(1+(block_k_idx-1)*(p+1)+1);
        
        [ H_bar_k_factor, P_k_vec, change_k ] = SYMUPD_mex(H_bar_k_factor, P_k_vec, change_k, b_k, u_tilde_k, p);
        [ H_bar_k_factor, P_k_vec, change_k ] = SYMUPD_mex(H_bar_k_factor, P_k_vec, change_k, -b_k, v_tilde_k, p);
        
        % restore L and B
        [L_k, B_k] = RESTORE(H_bar_k_factor, change_k, p);
        
        % update B and get B_sqrt from f(B, change)
        [B_k_sqrt_V, B_k_sqrt_D_vec] = get_B_sqrt_eig(B_k, change_k, eig_LB, p);
        
        %% record result
        theta_hat_ks_all(:, k+1, rep_idx) = theta_hat_k;
        loss_ks_all(1, k+1, rep_idx) = loss(theta_hat_k);
    end
end

theta_hat_ks = mean(theta_hat_ks_all, 3); % dim: p-by-n
loss_ks = mean(loss_ks_all, 3); % dim: 1-by-n
end

function d = get_delta(P_vec, L, V, D_vec, temp, p)
% get d such that (P' * L * V) * D * (V' * L' * P) * d = temp
P_d = L \ temp(P_vec); % V * D * V' * L' * (P * d) = inv(L) * (P * temp)
P_d = V' * P_d; % D * V' * L' * (P * d) = V' * inv(L) * (P * temp)
P_d = P_d ./ D_vec; % V' * L' * (P * d) = inv(D) * V' * inv(L) * (P * temp)
P_d = V * P_d; % M' * (P * d) = V * inv(D) * V' * inv(L) * (P * temp)
% Q_d = V * ((V' * P_d) ./ D_vec); % above three steps in one line

P_d = L' \ P_d; % (P * d) = inv(L') * V * inv(D) * V' * inv(L) * (P * temp)

Pt_vec = zeros(1,p); Pt_vec(P_vec) = 1:p; % vector form of P'
d = P_d(Pt_vec); % d = Q' * inv(L') * V * inv(D) * V' * inv(L) * (P * temp)
end

function td = get_tdelta(P_vec, L, V, D_vec, temp, p)
% get td such that (P' * L * V) * D * (V' * L' * P) * temp = td
td = V \ (L \ temp(P_vec)); % V' * L' * (P * temp)
td = D_vec .* td; % D * V' * L' * (P * temp)
td = L * (V * td); % L * V * D * V' * L' * (P * temp)
Pt_vec = zeros(1,p); Pt_vec(P_vec) = 1:p; % vector form of P'
td = td(Pt_vec); % P' * L * V * D * V' * L' * (P * temp)
end

function [V, D_vec] = get_B_sqrt_eig(B, change, eig_threshold, p)
V = zeros(p,p);
D_vec = zeros(p,1);

i = 1;
while i <= p
    block_size = change(i);
    if block_size == 1
        V(i,i) = 1;
        D_vec(i) = sqrt(max(eig_threshold, abs(B(i,i))));
        i = i + 1;
    else
        [block_V, block_D] = eig(B(i:i+1, i:i+1));
        D_vec(i:i+1) = sqrt(max(eig_threshold, abs(diag(block_D))));
        V(i:i+1,i:i+1) = block_V;
        i = i + 2;
    end
end
end


function [ M, D ] = RESTORE( A, CHANGE, p )
block_idx = find(CHANGE == 2);

M = tril(A); % get the lower triangular part
M(1:(p+1):(p^2)) = 1;
M(1+(block_idx-1)*(p+1)+1) = 0;

D = zeros(p); % block diagonal
D(1:(p+1):(p^2)) = A(1:(p+1):(p^2));
D(1+(block_idx-1)*(p+1)+1) = A(1+(block_idx-1)*(p+1)+1);
D(1+(block_idx-1)*(p+1)+p) = A(1+(block_idx-1)*(p+1)+1);

% for i = 1:p
%     if CHANGE(i) == 2
%         M(i+1,i) = 0;
%         M(i+1,i+1) = 1;
%         
%         D(i+1,i) = A(i+1,i);
%         D(i,i+1) = A(i+1,i);
%         D(i+1,i+1) = A(i+1,i+1);
%         
%         i = i + 1;
%     end
% end
end
