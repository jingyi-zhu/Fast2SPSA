function A = SCALE( A, change, scaler, p )
%A = scaler * A;
for i = 1:p
    A(i,i) = scaler * A(i,i);
    if (change(i) == 2)
        A(i+1,i) = scaler * A(i+1,i);
    end
end
