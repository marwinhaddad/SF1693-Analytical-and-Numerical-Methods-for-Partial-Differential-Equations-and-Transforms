function Xn = get_fixpoint(x)
    global h F tol
    diff = 1;
    
    while diff > tol
        Xn = x + h .* F(x);
        diff = abs(x - Xn);
        x = Xn;
    end
end
