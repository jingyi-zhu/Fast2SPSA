N = 100;
diff = zeros(N,1);
p = 0.01;
for n = 1 : 1 : N
    
    x = 0 : 1 : n;
    y_bin = binocdf(x,n,p);
    y_nor = normcdf(x, n*p, sqrt(n*p*(1-p)));
    
    % scatter(x, y_bin - y_nor);
    
    
    diff(n) = max(abs(y_bin - y_nor));
end
X = 1 : N; Y = 1 ./ X;
plot(X, diff, X, Y); 
legend('true','bound')
% max(diff) - 1/n


% p = 0.1;
% x = [0 : 1 : n];
% y_bin = binocdf(x,n,p);
% y_nor = normcdf(x, n*p, sqrt(n*p*(1-p))); 
% y_bin - y_nor
% diff = max(abs(y_bin - y_nor));