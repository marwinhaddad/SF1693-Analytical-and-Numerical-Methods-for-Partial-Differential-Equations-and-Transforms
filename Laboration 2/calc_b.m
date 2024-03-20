function b = calc_b(xs, fun, N, n)
    b = zeros(N, 1);
    
    for i = 1:N
        b(i) = quadr_b(xs(i), xs(i+1), n, fun, i) + quadr_b(xs(i+1), xs(i+2), n, fun, i);
    end
end