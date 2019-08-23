function gradient_fn = skewed_quartic_gradient(theta, B)
gradient_fn = B' * (2*B*theta + 0.3*sum((B*theta).^2) + 0.04*sum((B*theta).^3));
