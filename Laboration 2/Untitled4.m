
xs = linspace(0, 1, 10+1);
phi1 = phi(2, 0.05, xs)


function value = phi(j, x, xs)
        if x <= xs(j-1) || x >= xs(j+1)
            value = 0;
        elseif x >= xs(j-1) && x <= xs(j)
            value = (x-xs(j-1))./(xs(j)-xs(j-1));
        elseif x >= xs(j) && x <= xs(j+1)
            value = (xs(j+1)-x)./(xs(j+1)-xs(j));
        end
end