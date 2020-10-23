function [R, N] = myDecomp(A)
%MYDECOMP calculates the orthonormal basis of the range space of matrix A
%and the orthonormal basis of the null space of matrix A.
%   This function uses rref to get the reduced row echelon form of A or A',
%   which gives the basis of the null space and the range space
%   respectively. Then, this function uses Gram-Schimidt process to get the
%   orthonormalized result.
%
%   author Zhirui Dai
%   date October 22 2020

    [n, m] = size(A);
    
    % check dimension
    if n == 0 || m == 0
        R = zeros(n, 0);
        N = zeros(m, 0);
        disp('Zero dimension is not allowed!');
        return
    end
    
    % range space: {x in R^n: Aw = x for any w in R^m}
    [R, jb] = rref(A'); % jb is a vector containing the index of basis row
    R = R(jb, :)';      % basis row of A' is the basis column of A
    r = length(jb);
    [R, ~] = gram_schmidt(R);

    % null space: {x in R^m: Ax = 0}
    [Ntmp, jb] = rref(A);
    N = zeros(m, m - r);
    cnt = 1;
    % a method to quickly copy the numbers from Ntmp to N
    for i = 1 : m
        if ~any(jb == i)     % i is not in jb
            N(i, cnt) = 1;
            for j = 1 : length(jb)
                N(jb(j), cnt) = -Ntmp(j, i);
            end
            cnt = cnt + 1;
        end
    end
    N = gram_schmidt(N);
    
    % examine the answer
    if ~(rank([R orth(A)]) == r)
        disp('test of range space and orth failed');
    end
    
    if ~(rank([N null(A)]) == m - r)
        disp('test of nullspace and null failed');
    end
end
