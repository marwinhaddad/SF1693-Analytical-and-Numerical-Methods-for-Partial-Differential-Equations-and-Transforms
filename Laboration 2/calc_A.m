function A = calc_A(xs, fun, N, n)
    A = zeros(N, N);

    for i = 1:N
        for j = 1:N
            A(j, i) = quadr_A(xs(i), xs(i+1), n, fun, i, j) + quadr_A(xs(i+1), xs(i+2), n, fun, i, j);
        end
    end
end