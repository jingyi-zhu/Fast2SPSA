function loss_fn = skewed_quartic_loss(theta, B)
loss_fn = theta'*(B'*B)*theta + 0.1*sum((B*theta).^3) + 0.01*sum((B*theta).^4) + normrnd(0,0.05);
% normrnd(mu,sigma)