function E = get_energy(x)
    global N
    H = @(q1, q2, p1, p2) 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2);
    E = zeros(N+1, 1);
    for i = 1:N+1
        E(i) = H(x(i,1), x(i,2), x(i,3), x(i,4));
    end
end