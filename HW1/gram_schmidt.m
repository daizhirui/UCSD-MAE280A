function [Q, R] = gram_schmidt(A)
%GRAM_SCHIMIDT calculates the orthonormalization of columns in A, which
%gives the result of QR decomposition.
%   Suppose the size of A is n-by-m and its rank is r.
%   Q is a n-by-r matrix containing r orthonormalized column vectors. R is
%   a coefficient matrix of size r-by-m, where every column is the
%   projection of the corresponding colum of A onto Q.
%
%   author Zhirui Dai
%   date October 22 2020

    [n, m] = size(A);
    
    Q = zeros(n, m);
    R = zeros(m, m);
    
    if n == 0 || m == 0
        disp('Zero dimension is not allowed!');
        return
    end
    
    % must go through all columns because somebody may be evil and hides
    % those independent columns in the back of A
    for i = 1 : m
        R(1:i-1, i) = Q(:, 1:i-1).' * A(:, i);      % get projection
        v = A(:, i) - Q(:, 1:i-1) * R(1:i-1, i);    % orthogonalize
        R(i, i) = norm(v);                          % normalize
        if R(i, i) < 1.0E-16
            continue    % a dependent vector or a zero vector!
        end
        Q(:, i) = v / R(i, i);
    end
    
    % the Gram Schimidt process is done, but we can remove some zero
    % columns in Q
    indices = zeros(m);
    rank = 0;
    for i = 1 : m
        if any(abs(Q(:, i)) > 0)
            rank = rank + 1;
            indices(rank) = i;
        end
    end
    Q = Q(:, indices(1:rank));
    R = R(indices(1:rank), :);
end

