% late updated: 2019_01_24
% compare SL_2SP, FA_2SP

clear all;
p = 100;
seed = 100;

rep = 20;

% main parameters
n = 50000; % number of iterations
loss_ks = nan(4,n,rep);

% parameters for gain sequences
alpha = 0.602; gamma = 0.101;
% alpha = 1; gamma = 1/6;
a = 0.04; A = 1000; c = 0.05; c_tilde = c;
% a = 100; A = 100; c = 0.05; c_tilde = 0.05;
w = 0.01; d = 0.501;

theta_blocking_threshold = 1;

% generate gain, perturbation, weight sequences
% for k = 0
a_0 = a / (1+A)^alpha; c_0 = c; c_tilde_0 = c_tilde; w_0 = 1;
% for k = 1,...,n
a_ks = a ./ (((1:n) + 1 + A) .^ alpha); % a_k = a/(k+1+A)^alpha
c_ks = c ./ (((1:n) + 1) .^ gamma); % c_k = c/(k+1)^gamma
c_tilde_ks = c_tilde ./ (((1:n) + 1) .^ gamma); % c_tilde_k = c_tilde/(k+1)^gamma
% w_ks = 1 ./ ((1:n) + 1); % w_k = 1/(k+1)

% loss function and optimal values
B_quartic = triu(ones(p)) / p; % the matrix used in skewed quartic function
theta_star = zeros(p, 1);
loss_star = skewed_quartic_loss_noise_free(theta_star, B_quartic);

% parameters space for theta
bound = 10;
theta_min = - bound * ones(p,1); % lower bounds on theta
theta_max = bound * ones(p,1); % upper bounds on theta

% initialization
theta_hat_0 = 1 * ones(p,1);
loss_0 = skewed_quartic_loss_noise_free(theta_hat_0, B_quartic);
H_bar_0 = eye(p);

% Q * H * Q' = M * B * M'
% Q: permutation matrix
% H: estimated Hession
% M: lower triangle matrix
% B: block diagonal matrix of size 2
Q_0_vec = 1:p; % store Q in vector form
change_0 = ones(1,p); % block structure indicator for B
H_bar_0_factor = eye(p); % store M and B

for algo_idx = [1,2] % SL_2SP, FA_2SP, SL_E2SP, FA_E2SP
    rng(seed);
    
    for rep_idx = 1:rep
        %% initialization: case k == 0
        delta_0 = 2 * round(rand(p, 1)) - 1;
        % gradient estimates (two-sided)
        theta_hat_0_plus = theta_hat_0 + c_0 * delta_0;
        theta_hat_0_minus = theta_hat_0 - c_0 * delta_0;
        % estimate gradient if using 2SPSA
        y_hat_0_plus = skewed_quartic_loss(theta_hat_0_plus, B_quartic);
        y_hat_0_minus = skewed_quartic_loss(theta_hat_0_minus, B_quartic);
        G_hat_0 = (y_hat_0_plus - y_hat_0_minus) / (2 * c_0) * delta_0;
        % Hessian estimate
        H_bar_k = H_bar_0;
        if (algo_idx == 2) || (algo_idx == 4) % for FA_2SP and FA_E2SP
            Q_k_vec = Q_0_vec;
            H_bar_k_factor = H_bar_0_factor;
            change_k = change_0;
        end
        s_0 = H_bar_0 \ (-G_hat_0);

        %% update theta_hat_k
        % blocking
        theta_hat_new = theta_hat_0 + a_0 * s_0; % theta_hat_{k=1}
        if norm(theta_hat_new - theta_hat_0) < theta_blocking_threshold
            theta_hat_k = theta_hat_new;
            theta_hat_k = min(theta_hat_k, theta_max);
            theta_hat_k = max(theta_hat_k, theta_min);
        else
            theta_hat_k = theta_hat_0;
        end
        loss_ks(algo_idx, 1, rep_idx) = (skewed_quartic_loss_noise_free(theta_hat_k, B_quartic) - loss_star) / (loss_0 - loss_star);

        %% iteration starts
        for iter_idx = 1:(n-1)
            % gain sequences
            a_k = a_ks(iter_idx);
            c_k = c_ks(iter_idx);
            c_tilde_k = c_tilde_ks(iter_idx);

            % perturbations
            delta_k = 2 * round(rand(p, 1)) - 1;
            delta_tilde_k = 2 * round(rand(p, 1)) - 1;

            % adaptive weight
            if (algo_idx == 1) || (algo_idx == 2)
                w_k = w / (iter_idx + 1).^d;
            elseif (algo_idx == 3) || (algo_idx == 4)
                % optimal weight
    %             w_k = c_k^2 * c_tilde_k^2 / ...
    %                 (c_0^2*c_tilde_0^2 + sum(c_ks(1:iter_idx).^2 .* c_tilde_ks(1:iter_idx).^2));
                w_k = w / (iter_idx + 1).^d;
            end

            %% gradient estimates (two-sided)
            % generate pertubations
            theta_hat_k_plus = theta_hat_k + c_k * delta_k;
            theta_hat_k_minus = theta_hat_k - c_k * delta_k;

            % estimate gradient if using 2SP
            y_hat_k_plus = skewed_quartic_loss(theta_hat_k_plus, B_quartic);
            y_hat_k_minus = skewed_quartic_loss(theta_hat_k_minus, B_quartic);
            G_hat_k = (y_hat_k_plus - y_hat_k_minus) / (2 * c_k) * delta_k;

            %% Hessian estimates (one-sided)
            theta_tilde_k_plus = theta_hat_k_plus + c_tilde_k * delta_tilde_k;
            theta_tilde_k_minus = theta_hat_k_minus + c_tilde_k * delta_tilde_k;
            y_tilde_k_plus = skewed_quartic_loss(theta_tilde_k_plus, B_quartic);
            y_tilde_k_minus = skewed_quartic_loss(theta_tilde_k_minus, B_quartic);
            if (algo_idx == 1) || (algo_idx == 2) % SL_2SP or FA_2SP
                d_k = 1 - w_k;
                delta_y_k = (y_tilde_k_plus - y_hat_k_plus) - (y_tilde_k_minus - y_hat_k_minus);
                b_k = w_k * delta_y_k / (4 * c_k * c_tilde_k);
                u_k = delta_tilde_k; v_k = delta_k;
                u_k_norm = norm(u_k); v_k_norm = norm(v_k);
                u_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k + u_k_norm/v_k_norm*v_k);
                v_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k - u_k_norm/v_k_norm*v_k);            
            elseif (algo_idx == 3) || (algo_idx == 4) % SL_E2SP or FA_E2SP
                d_k = 1;
                delta_y_k = (y_tilde_k_plus - y_hat_k_plus) - (y_tilde_k_minus - y_hat_k_minus);
                b_k = w_k * (delta_y_k / (2 * c_k * c_tilde_k) - delta_k'* H_bar_k * delta_tilde_k) / 2;
                u_k = delta_tilde_k; v_k = delta_k;
                u_k_norm = norm(u_k); v_k_norm = norm(v_k);
                u_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k + u_k_norm/v_k_norm*v_k);
                v_tilde_k = sqrt(v_k_norm/(2*u_k_norm)) * (u_k - u_k_norm/v_k_norm*v_k);
            end

            %% quasi-Newton step
            if (algo_idx == 1) || (algo_idx == 3) % SL_2SP or SL_E2SP
                H_bar_k = d_k * H_bar_k + b_k * (u_tilde_k * u_tilde_k' - v_tilde_k * v_tilde_k');
            elseif (algo_idx == 2) || (algo_idx == 4) % FA_2SP or FA_E2SP
                % Hbar = d_k * H_bar;
                H_bar_k_factor(1:(p+1):(p^2)) = d_k * H_bar_k_factor(1:(p+1):(p^2));
                block_k_idx = find(change_k == 2);
                H_bar_k_factor(1+(block_k_idx-1)*(p+1)+1) = d_k * H_bar_k_factor(1+(block_k_idx-1)*(p+1)+1);
                % Hbar = Hbar + b_k * (uktilde * uktilde')
                [H_bar_k_factor, Q_k_vec, change_k] = SYMUPD(H_bar_k_factor, Q_k_vec, change_k, b_k, u_tilde_k, p);
                % Hbar = Hbar - b_k * (vktilde * vktilde')
                [H_bar_k_factor, Q_k_vec, change_k] = SYMUPD(H_bar_k_factor, Q_k_vec, change_k, -b_k, v_tilde_k, p);
                if algo_idx == 4 % FA_E2SP
                    H_bar_k = d_k * H_bar_k + b_k * (u_tilde_k * u_tilde_k' - v_tilde_k * v_tilde_k');
                end
            end

            %% modified-Newton step
            if (algo_idx == 1) || (algo_idx == 3) % SL_2SP or SL_E2SP
                % make it positive definite
                H_barbar_k = sqrtm(H_bar_k * H_bar_k' + 10^(-4) * exp(-iter_idx) * eye(p));

            elseif (algo_idx == 2) || (algo_idx == 4) % FA_2SP or FA_E2SP
                % get the factorization matrix from H_bar_k_factor
                [M_k, B_k] = RESTORE(H_bar_k_factor, change_k, p);
                [V_k, Lambda_k] = eig(B_k); % B = V * Lambda * V'
                eig_vec = diag(Lambda_k);
    %             eig_bar_vec = max(eig_vec, 10^(-4));
                eig_bar_vec = max(max(abs(eig_vec), 10^(-8)),  10^(-8) * p * max(abs(eig_vec)));
            end
            %% descent direction
            if (algo_idx == 1) || (algo_idx == 3) % SL_2SP or SL_E2SP
                s_k = H_barbar_k \ (-G_hat_k);
            elseif (algo_idx == 2) || (algo_idx == 4) % FA_2SP or FA_E2SP
                % descent direction: (Q' * M * V) * D * (V' * M' * Q) * s = -g
                % M * V * D * V' * M' * (Q*s) = -(Q*g)
                Qs_k = M_k \ -G_hat_k(Q_k_vec); % V * Lambda * V' * M' * (Q*s) = inv(M) * (-(Q*g))
                Qs_k = V_k' * Qs_k; % D * V' * M' * (Q*s) = V' * [inv(M) * (-(Q*g))]
                Qs_k = Qs_k ./ eig_bar_vec; % V' * M' * (Q*s) = inv(Lambda) * [V' * inv(M) * (-(Q*g))]
                Qs_k = V_k * Qs_k; % M' * (Q*s) = V * [inv(Lambda) * V' * inv(M) * (-(Q*g))]
                Qs_k = M_k' \ Qs_k; % (Q*s) = inv(M') * [V * inv(Lambda) * V' * inv(M) * (-(Q*g))]
                trans_Q_k_vec = zeros(1,p); trans_Q_k_vec(Q_k_vec) = 1:p; % vector form of Q'
                s_k = Qs_k(trans_Q_k_vec); % s = Q' * [inv(M') * V * inv(Lambda) * V' * inv(M) * (-(Q*g))]
            end

            %% update theta_hat_k
            % blocking
            theta_hat_k_new = theta_hat_k + a_k * s_k;
            if norm(theta_hat_k_new - theta_hat_k) < theta_blocking_threshold
                theta_hat_k = theta_hat_k_new;
                theta_hat_k = min(theta_hat_k, theta_max);
                theta_hat_k = max(theta_hat_k, theta_min);
            end
            loss_ks(algo_idx, iter_idx+1, rep_idx) = (skewed_quartic_loss_noise_free(theta_hat_k, B_quartic) - loss_star) / (loss_0 - loss_star);        
        end
    end
end

loss_ks_avg_rep = mean(loss_ks, 3);

figure;
hold on
plot(0:n, [1, loss_ks_avg_rep(1, :)], 'black-', 'Linewidth', 1.5);
plot(0:n, [1, loss_ks_avg_rep(2, :)], 'black--', 'Linewidth', 1.5);
% plot(0:n, [1, loss_ks_avg_rep(3, :)], 'black:', 'Linewidth', 1.5);
% plot(0:n, [1, loss_ks_avg_rep(4, :)], 'black-.', 'Linewidth', 1.5);

title({'Normalized Loss Function', 'Value Comparison'}, 'FontSize', 18)
% set(gca,'FontSize',20,'YScale','log')
set(gca,'FontSize',20)
xlabel({'Iteration $k$'}, 'Interpreter', 'latex', 'FontSize', 18)
ylabel({'Average Normalized Loss', 'Function Value per Replicate'}, 'FontSize', 18)
% legend({'Original 2SPSA','Efficient 2SPSA','Original E2SPSA','Efficient E2SPSA'},'FontSize',13)
legend({'Original 2SPSA','Efficient 2SPSA'},'FontSize',13)
hold off

save('compare_2SPSA_loss_per_iteration_2019_5_27');

fig = gcf;
% set(gcf, 'InvertHardCopy', 'off');
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'compare_2SPSA_loss_per_iteration_2019_5_27','-dpdf');