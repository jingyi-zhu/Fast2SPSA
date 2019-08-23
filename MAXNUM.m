function [ MAXMUM, I ] = MAXNUM( A, B, C )
% MAXABSINDX
I = 1;
T = abs(A);
S = abs(B);
if (S <= T)
    % Go to 10
else
    T = S;
    I = 2;
end
% 10
S = abs(C);
if (S <= T)
    % Go to 20
else
    T = S;
    I = 3;
end
% 20
MAXMUM = T;
end

