function M = complex2cols(C)
    a = real(C);
    b = imag(C);
    ncol = size(a, 2);
    for k = 1:ncol
        M(:,2*k-[1 0]) = [a(:,k) b(:,k)];
    end
