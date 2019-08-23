function hessian_fn = skewed_quartic_hessian(theta, B)
hessian_fn = B' * diag(2 + 0.6*B*theta + 0.12*sum((B*theta).^2)) * B;