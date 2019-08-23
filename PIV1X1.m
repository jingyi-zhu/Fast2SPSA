function [ D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, Q1, Q2, Q3, I1, K, N, KP1, KP2 ]...
    = PIV1X1(D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, Q1, Q2, Q3, I1, K, N)
%PIV1X1 Summary of this function goes here
%   Detailed explanation goes here
KP1 = K + 1;
KP2 = K + 2;
PIV1X1CASENUM = I1 * 10;
while true
    switch (PIV1X1CASENUM)
        case 20 % The max element is D22
            T = D11;
            D11 = D22;
            D22 = T;
            
            T = D32;
            D32 = D31;
            D31 = T;
            
            T = B2;
            B2 = B1;
            B1 = T;
            Q1 = 2;
            Q2 = 1;
            
            % Go to 10
            PIV1X1CASENUM = 10;
            continue;
            
        case 30 % The max element is D33
            T = D11;
            D11 = D33;
            D33 = T;
            
            T = D21;
            D21 = D32;
            D32 = T;
            
            T = B1;
            B1 = B3;
            B3 = T;
            
            Q1 = 3;
            Q3 = 1;
            % Go to 10
            PIV1X1CASENUM = 10;
            continue;
            
        case 10 % The max element is D11
            D22 = D22 - (D21 * D21) / D11;
            D32 = D32 - (D31 * D21) / D11;
            D33 = D33 - (D31 * D31) / D11;
            B1 = B1 / D11;
            B2 = B2 - B1 * D21;
            B3 = B3 - B1 * D31;
            D21 = D21 / D11;
            D31 = D31 / D11;
            CHANGE(K) = 1;
            break;
    end
end
end
