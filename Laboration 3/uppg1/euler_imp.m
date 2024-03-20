function x = euler_imp()
    global N h X F
    x = X;
    
    for i = 1:N
        x(i+1, :) = get_fixpoint(x(i, :));
        x(i+1, :) = x(i, :) + h .* F(x(i+1, :));
    end
end