function value = quadr_A(a, b ,n, fun, i, j)
    if n == 1
        value = (b-a)*fun((b+a)/2, i, j);
    else    
        value = ((b-a)/2)*fun(((b-a)/2)*(-1/sqrt(3))+(b+a)/2, i, j) + ((b-a)/2)*fun(((b-a)/2)*(1/sqrt(3))+(b+a)/2, i, j);
    end
end