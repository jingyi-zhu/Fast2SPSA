function loss_fn = skewed_quartic_loss_noise_free(theta, B)
loss_fn = theta'*(B'*B)*theta + 0.1*sum((B*theta).^3) + 0.01*sum((B*theta).^4);
% normrnd(mu,sigma)