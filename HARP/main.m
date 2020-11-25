clear all; 
rng(19);
algo_list = [1 2];

%% initialization
p = 20; loss = @(theta) skewed_quartic(theta,p); theta_star = zeros(p, 1); 
 
range = 20;
theta_0 = rand(p, 1) * 2 * range - range; 
loss_0 = loss(theta_0); loss_star = loss(theta_star);

%% noise
noise_var = 0.1;
loss_noisy = @(theta) skewed_quartic(theta,p) + mvnrnd(zeros(1,1),  noise_var *  eye(1) )';

%% itration
n = 2e4; rep = 25;  alpha = 0.602; gamma = 0.101;
A = 200;

%% 2nd-order SPSA
a_SPSA_second = 7;  
c_SPSA_second = 1;
if ismember(1, algo_list)
    [SPSA_second_theta_hat_ks, SPSA_second_loss_ks] = ...
        SPSA_second(a_SPSA_second,A,alpha,c_SPSA_second,gamma,loss,loss_noisy,n,rep,theta_0);
    SPSA_second_theta_norm = norm(SPSA_second_theta_hat_ks(:, end))
else
    SPSA_second_loss_ks = nan(1,n);
    SPSA_second_theta_hat_ks = nan(p,n);
end
SPSA_second_loss_ks_normalized = (SPSA_second_loss_ks - loss_star) / (loss_0 - loss_star);
SPSA_second_theta_ks_normalized = vecnorm(SPSA_second_theta_hat_ks - theta_star * ones(1,n)) / norm(theta_0 - theta_star);


%% 2nd-order HARP
a_HARP_second = 10;  
c_HARP_second = 1;
if ismember(2, algo_list)
    [HARP_second_theta_hat_ks, HARP_second_loss_ks] = ...
        HARP_second(a_HARP_second,A,alpha,c_HARP_second,gamma,loss,loss_noisy,n,rep,theta_0);
    HARP_second_theta_norm = norm(HARP_second_theta_hat_ks(:, end))
else
    HARP_second_loss_ks = nan(1,n);
    HARP_second_theta_hat_ks = nan(p,n);
end
HARP_second_loss_ks_normalized = (HARP_second_loss_ks - loss_star) / (loss_0 - loss_star);
HARP_second_theta_ks_normalized = vecnorm(HARP_second_theta_hat_ks - theta_star * ones(1,n)) / norm(theta_0 - theta_star);

%% plot
figure(1); clf;
hold on
plot(0:n, ([1, SPSA_second_loss_ks_normalized]), '--r', 'Linewidth', 1.5);
plot(0:n, ([1, HARP_second_loss_ks_normalized]), '-b', 'Linewidth', 1.5);

set(gca,'FontSize',20, 'YScale', 'log')
legend({'SPSA', 'HARP'},'FontSize',13)
title({'Test Function: Skew-Quartic',  ['noise-var =', num2str(noise_var), ', dim = ', num2str(p), ', init-range: [-', num2str(range), ',', num2str(range), '],  rep = ', num2str(rep)], ['a-SPSA = ',  num2str(a_SPSA_second), ', c-SPSA=',  num2str(c_SPSA_second) ], ['a-HARP= ',  num2str(a_HARP_second), ', c-HARP=',  num2str(c_HARP_second)] }, 'FontSize', 18)
xlabel({'Iteration $k$'}, 'Interpreter', 'latex', 'FontSize', 15)
ylabel({'Average Normalized Loss', 'Function Value per Replicate'}, 'FontSize', 15)
hold off

figure(2); clf;
hold on
plot(0:n, ([1, SPSA_second_theta_ks_normalized]), '--r', 'Linewidth', 1.5);
plot(0:n, ([1, HARP_second_theta_ks_normalized]), '-b', 'Linewidth', 1.5);

set(gca,'FontSize',20, 'YScale', 'log')
legend({'SPSA', 'HARP'},'FontSize',13)
title({'Test Function: Skew-Quartic',  ['noise-var =', num2str(noise_var), ', dim = ', num2str(p), ', init-range: [-', num2str(range), ',', num2str(range), '],  rep = ', num2str(rep)], ['a-SPSA = ',  num2str(a_SPSA_second), ', c-SPSA=',  num2str(c_SPSA_second) ], ['a-HARP= ',  num2str(a_HARP_second), ', c-HARP=',  num2str(c_HARP_second)] }, 'FontSize', 18)
xlabel({'Iteration $k$'}, 'Interpreter', 'latex', 'FontSize', 15)
ylabel({'Average Normalized Distance to Optimum'}, 'FontSize', 15)
hold off



function loss_fn = skewed_quartic(theta,p)
B = triu(ones(p)) / p;
loss_fn = theta'*(B'*B)*theta + 0.1*sum((B*theta).^3) + ...
    0.01*sum((B*theta).^4);
end

function loss_fn = Rosenbrock(theta,p)
loss_fn = 0;
for i = 1:(p-1)
    loss_fn = loss_fn + 100 * (theta(i+1)^2 - theta(i))^2 + ...
        (theta(i) - 1)^2;
end
end
