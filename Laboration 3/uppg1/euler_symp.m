function x = euler_symp()
    global N h X f
    Q = X(:, 1:2);
    P = X(:, 3:4);
    
    for i = 1:N
        Q(i+1, :) = Q(i, :) + h .* P(i, :);
        P(i+1, :) = P(i, :) + h .* f(Q(i+1, :));
    end
    x = [Q P];
end