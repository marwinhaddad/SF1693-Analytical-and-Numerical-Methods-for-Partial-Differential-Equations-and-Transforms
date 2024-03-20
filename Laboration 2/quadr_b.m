function value = quadr_b(a, b, n, fun, i)
    if n == 1
        value = (b-a)*fun((b+a)/2, i);
    else    
        value = ((b-a)/2)*fun(((b-a)/2)*(-1/sqrt(3))+(b+a)/2, i) + ((b-a)/2)*fun(((b-a)/2)*(1/sqrt(3))+(b+a)/2, i);
    end
end