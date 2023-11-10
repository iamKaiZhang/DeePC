function H = make_hankel(x, L)
    [m, K] = size(x);
    H = zeros(m*L, K-L+1);

    for i = 1:K-L+1
        H(:, i) = reshape(x(:, i:i+L-1), [], 1);
    end
end