function x = euler_exp()
    global N h X F
    x = X;

    for i = 1:N
        x(i+1, :) = x(i, :) + h.* F(x(i, :));
    end
end