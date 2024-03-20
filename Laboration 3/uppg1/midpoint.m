function x = midpoint()
    global N h X F f
    x = X;
    
    for i = 1:N
        x(i+1, :) = get_fixpoint(x(i, :));
        x(i+1, :) = x(i, :) + h .* F(1/2 .* (x(i, :) + x(i+1, :)));
    end
end