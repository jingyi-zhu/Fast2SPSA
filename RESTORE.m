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
