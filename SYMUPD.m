function [A, Q, CHANGE] = SYMUPD(A, Q, CHANGE, SIGMA, Z, p)
ALEA = (1 + sqrt(17))/8;
%ALEA = 0;
W = NaN(p, 1);
W(1:p) = Z(Q(1:p)); % W = Q * Z
K = 1;

D11 = 0; D21 = 0; D22 = 0;
B1 = 0; B2 = 0;
KP1 = 0; KP2 = 0;

CASENUM = 101;
while true
    switch(CASENUM) % 101, 103, 108, 115, 117
        case 101
            % The processing begins here
            KP1 = K + 1;
            KP2 = K + 2;
            if (K < p) && (CHANGE(KP1) < 0)
                % Go to 115
                CASENUM = 115;
                continue;
            else
                % The next block is a 1X1 block
                T = W(K);
                B1 = SIGMA * T;
                D11 = A(K,K) + B1 * T;
                if KP1 > p
                    % Go to 202
                else
                    % Do 102
                    % J = KP1:p;
                    W(KP1:p) = W(KP1:p) - A(KP1:p,K) * T;
                end
                % 202
                % Go to 103
                CASENUM = 103;
                continue;
            end
        case 103
            % Enter 1X1
            if K < p
                % Go to 104
            else
                A(K,K) = D11;
                % The decomposition is complete if K = N
                break;
            end
            % 104
            U1 = abs(D11);
            U0 = abs(B1);
            
            if (U1 < (ALEA*U0)) && ((U1*abs(SIGMA)) <= (ALEA*U0*U0))
                % Go to 106
            else
                SIGMA = SIGMA - B1 * B1 / D11;
                B1 = B1 / D11;
                A(K,K) = D11;
                % Update the K-th column of M
                if KP1 > p
                    % Go to 205
                else
                    % Do 105
                    % J = KP1:p;
                    A(KP1:p,K) = A(KP1:p,K) + B1 * W(KP1:p);
                end
                % 205
                K = KP1;
                % Go to 101
                CASENUM = 101;
                continue;
            end
            % 106
            if (KP2 <= p) && (CHANGE(KP2) <= 0)
                % Go to 108
                CASENUM = 108;
                continue;
            else
                % A 2X2 block is formed by combining the next 1X1 block with block K
                B2 = W(KP1);
                T = B2;
                D21 = B2 * B1;
                B2 = SIGMA * B2;
                D22 = A(KP1,KP1) + B2 * T;
                L1 = A(KP1,K);
                D22 = D22 + L1 * D21;
                B2 = B2 + L1 * B1;
                D21 = D21 + L1 * D11;
                D22 = D22 + L1 * D21;
                
                % Include information from the (K+1)-st column of M
                if KP2 > p
                    % Go to 207
                else
                    % Do 107
                    % J = KP2:p;
                    W(KP2:p) = W(KP2:p) - A(KP2:p,KP1) * T;
                    A(KP2:p,K) = A(KP2:p,K) - A(KP2:p,KP1) * L1;
                end
                % 207
                % Go to 117
                CASENUM = 117;
                continue;
            end
        case 108 % If this portion of code is reached we are in the case of a 1X1 singular block followed by a 2X2 block
            T1 = W(KP1);
            T2 = W(KP2);
            B2 = SIGMA * T1;
            B3 = SIGMA * T2;
            D22 = A(KP1,KP1) + B2 * T1;
            D32 = A(KP2,KP1) + B3 * T1;
            D33 = A(KP2,KP2) + B3 * T2;
            D21 = T1 * B1;
            D31 = T2 * B1;
            L1 = A(KP1,K);
            L2 = A(KP2,K);
            T = L2 * D11;
            D33 = D33 + L2 * (2*D31 + T);
            D31 = D31 + T;
            D32 = D32 + L1 * D31 + L2 * D21;
            T = L1 * D11;
            D22 = D22 + L1 * (2*D21 + T);
            D21 = D21 + T;
            B2 = B2 + L1 * B1;
            B3 = B3 + L2 * B1;
            KP3 = K + 3;
            % Include information from the (K+1)-st and (K+2)-nd columns of M
            if (KP3 > p) % Go to 209
            else % 109
                W(KP3:p) = W(KP3:p) - (A(KP3:p,KP1) * T1 + A(KP3:p,KP2) * T2);
                A(KP3:p,K) = A(KP3:p,K) - (A(KP3:p,KP1) * L1 + A(KP3:p,KP2) * L2);
            end
            % 209
            [U1, I1] = MAXNUM(D11, D22, D33);
            [U0, I0] = MAXNUM(D21, D31, D32);
            if (U1 < (ALEA * U0)) % Go to 112
            else % A 1X1 pivot will be used
                Q1 = 1; Q2 = 2; Q3 = 3;
                [D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, Q1, Q2, Q3, I1, K, p, KP1, KP2] ...
                    = PIV1X1(D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, Q1, Q2, Q3, I1, K, p);
                K1 = K - 1 + Q1;
                SIGMA = SIGMA - D11 * B1 * B1;
                % Update the K-th column of M
                if (KP3 > p) % Go to 210
                else
                    % 110
                    % J = KP3:p;
                    T = A(KP3:p,K);
                    A(KP3:p,K) = A(KP3:p,K1);
                    A(KP3:p,K1) = T;
                    A(KP3:p,K) = A(KP3:p,K) + D21 * A(KP3:p,KP1) + D31 * A(KP3:p,KP2) + B1 * W(KP3:p);
                end
                % 210
                KM1 = K - 1;
                % Intercahnge the corresponding rows of M
                if (KM1 < 1) % Go to 211
                else
                    % 111
                    % J = 1:KM1;
                    T = A(K,1:KM1);
                    A(K,1:KM1) = A(K1,1:KM1);
                    A(K1,1:KM1) = T;
                end
                % 211
                I = Q(K);
                Q(K) = Q(K1);
                Q(K1) = I;
                
                A(K,K) = D11;
                A(KP1,K) = D21;
                A(KP2,K) = D31;
                D11 = D22;
                D22 = D33;
                D21 = D32;
                B1 = B2;
                B2 = B3;
                K = KP1;
                KP1 = KP2;
                KP2 = KP2 + 1;
                % Go to 117
                CASENUM = 117;
                continue;
            end
            % 112, A 2X2 pivot will be used
            Q1 = 1;
            Q2 = 2;
            Q3 = 3; 
            [D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, SIGMA, Q1, Q2, Q3, I0, K, p, KP1, KP2] ...
                = PIV2X2(D11, D21, D31, D22, D32, D33, B1, B2, B3, CHANGE, SIGMA, Q1, Q2, Q3, I0, K, p);
            K1 = K - 1 + Q1;
            K2 = K - 1 + Q2;
            
            I = Q(K);
            Q(K) = Q(K1);
            Q(K1) = I;
            
            I = Q(KP1);
            Q(KP1) = Q(K2);
            Q(K2) = I;
            
            % Update the K-th and (K+1)-st colums of M
            if (KP3 > p) % Go to 213
            else
                % 113
                % J = KP3:p;
                T = A(KP3:p,K);
                A(KP3:p,K) = A(KP3:p,K1);
                A(KP3:p,K1) = T;
                
                T = A(KP3:p,KP1);
                A(KP3:p,KP1) = A(KP3:p,K2);
                A(KP3:p,K2) = T;
                
                A(KP3:p,K) = A(KP3:p,K) + D31 * A(KP3:p,KP2) + B1 * W(KP3:p);
                A(KP3:p,KP1) = A(KP3:p,KP1) + D32 * A(KP3:p,KP2) + B2 * W(KP3:p);
            end
            % 213, Interchange the correspondig rows of M
            KM1 = K - 1;
            % J = 1:KM1;
            T = A(K,1:KM1);
            A(K,1:KM1) = A(K1,1:KM1);
            A(K1,1:KM1) = T;

            T = A(KP1,1:KM1);
            A(KP1,1:KM1) = A(K2,1:KM1);
            A(K2,1:KM1) = T;
            
            % 114
            A(K,K) = D11;
            A(KP1,K) = D21;
            A(KP1,KP1) = D22;
            A(KP2,K) = D31;
            A(KP2,KP1) = D32;
            
            D11 = D33;
            B1 = B3;
            K = KP2;
            KP1 = K + 1;
            KP2 = K + 2;
            
            % Go to 103
            CASENUM = 103;
            continue;
        case 115 % The diagonal block beginning at entry K is 2X2
            T1 = W(K);
            T2 = W(KP1);
            B1 = SIGMA * W(K);
            B2 = SIGMA * W(KP1);
            D11 = A(K,K) + B1 * W(K);
            D21 = A(KP1,K) + B2 * W(K);
            D22 = A(KP1,KP1) + B2 * W(KP1);
            if (KP1 >= p) % Go to 117
            else
                % 116
                % J = KP2:p;
                W(KP2:p) = W(KP2:p) - (A(KP2:p,K) * T1 + A(KP2:p,KP1) * T2);
            end
            % 117
            CASENUM = 117;
            continue;
        case 117 % Enter 2X2
            T1 = 0;
            [U1, I1] = MAXNUM(D11, D22, T1);
            if (U1 >= (ALEA * abs(D21))) % Go to 119
            else
                % A 2X2 pivot will be used
                DET = D11 * D22 - D21 * D21;
                CHANGE(K) = 2;
                CHANGE(KP1) = DET;
                A(K,K) = D11;
                A(KP1,K) = D21;
                A(KP1,KP1) = D22;
                if (KP1 == p)
                    break;
                end
                T1 = (D22*B1 - D21*B2) / DET;
                T2 = (-D21*B1 + D11*B2) / DET;
                if (KP2 > p) % Go to 218
                else
                    % 118
                    % J = KP2:p;
                    T = W(KP2:p);
                    A(KP2:p,K) = A(KP2:p,K) + T1 * T;
                    A(KP2:p,KP1) = A(KP2:p,KP1) + T2 * T;
                end
                % 218
                SIGMA = SIGMA - (T1*B1 + T2*B2);
                K = KP2;
                
                % Go to 101
                CASENUM = 101;
                continue;
            end
            % 119 A 1X1 pivot will be used
            if I1 ~= 2
                % Go to 122
            else
                % Interchange the rows and columns of M if necessary
                T = D11;
                D11 = D22;
                D22 = T;

                T = B1;
                B1 = B2;
                B2 = T;

                I = Q(K);
                Q(K) = Q(KP1);
                Q(KP1) = I;

                if (KP2 > p)
                    % Go to 220
                else
                    % Do 120
                    % J = KP2:p;
                    T = A(KP2:p,K);
                    A(KP2:p,K) = A(KP2:p,KP1);
                    A(KP2:p,KP1) = T;
                end
                % 220
                KM1 = K - 1;
                if KM1 < 1
                    % Go to 221
                else
                    % Do 121
                    % J = 1:KM1;
                    T = A(K,1:KM1);
                    A(K,1:KM1) = A(KP1,1:KM1);
                    A(KP1,1:KM1) = T;
                end
                % 221
            end
            % 122 Process the two 1X1 blocks
            CHANGE(K) = 1;
            CHANGE(KP1) = 1;
            D22 = D22 - (D21*D21) / D11;
            D21 = D21 / D11;
            B2 = B2 - B1 * D21;
            B1 = B1 / D11;
            if KP2 > p
                % Go to 223
            else
                % Update the K-th column of M
                % Do 123
                % J = KP2:p;
                A(KP2:p,K) = A(KP2:p,K) + D21 * A(KP2:p,KP1) + B1 * W(KP2:p);
            end
            % 223
            A(K,K) = D11;
            A(KP1,K) = D21;
            SIGMA = SIGMA - B1 * D11 * B1;
            D11 = D22;
            B1 = B2;
            K = KP1;
            KP1 = KP2;
            KP2 = K + 2;
            
            % Go to 103
            CASENUM = 103;
            continue;
    end
end
end