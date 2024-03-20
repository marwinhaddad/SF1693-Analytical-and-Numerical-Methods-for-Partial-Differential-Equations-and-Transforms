function value = phi(i, x, xs)
    if i+2 > length(xs)
        i = i-1;
    end
    if x <= xs(i) || x >= xs(i+2)
        value = 0;
    elseif x >= xs(i) && x <= xs(i+1)
        value = (x-xs(i))./(xs(i+1)-xs(i));
    elseif x >= xs(i+1) && x <= xs(i+2)
        value = (xs(i+2)-x)./(xs(i+2)-xs(i+1));
    end
end