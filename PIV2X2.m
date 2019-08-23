function [ D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, SIGMA, Q1, Q2, Q3, I0, K, N, KP1, KP2 ]...
    = PIV2X2(D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, SIGMA, Q1, Q2, Q3, I0, K, N)
KP1 = K + 1;
KP2 = K + 2;
PIV2X2CASENUM = I0 * 10;
while true
    switch (PIV2X2CASENUM)
        case 20 % The max element is D31
            T = D22;
            D22 = D33;
            D33 = T;
            
            T = D21;
            D21 = D31;
            D31 = T;
            
            T = B2;
            B2 = B3;
            B3 = T;
            
            Q2 = 3;
            Q3 = 2;
            
            % Go to 10
            PIV2X2CASENUM = 10;
            continue;
            
        case 30 % The max element is D32
            T = D11;
            D11 = D22;
            D22 = D33;
            D33 = T;
            
            T = D21;
            D21 = D32;
            D32 = D31;
            D31 = T;
            
            T = B1;
            B1 = B2;
            B2 = B3;
            B3 = T;
            
            Q1 = 2;
            Q2 = 3;
            Q3 = 1;
            
            % Go to 10
            PIV2X2CASENUM = 10;
            continue;
            
        case 10 % The max element is D21, the 2X2 pivot is done here
            DET = D11 * D22 - D21 * D21;
            T = (D22*D31 - D21*D32) / DET;
            S = (-D21*D31 + D11*D32) / DET;
            B3 = B3 - (T*B1 + S*B2);
            D33 = D33 - (T*D31 + S*D32);
            D31 = T;
            D32 = S;
            T = (D22*B1 - D21*B2) / DET;
            S = (-D21*B1 + D11*B2) / DET;
            SIGMA = SIGMA - (T*B1 + S*B2);
            B1 = T;
            B2 = S;
            CHANGE(K) = 2;
            CHANGE(KP1) = DET;
            CHANGE(KP2) = 1;
            break;
    end
end
end