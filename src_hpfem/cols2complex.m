function M = cols2complex(A)
    ncol = size(A, 2)/2;
    for k = 1:ncol
        a = A(:, 2*k-1);
        b = A(:, 2*k);
        M(:, k) = complex(a, b);
    end