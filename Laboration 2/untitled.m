clear all;
clc;
close all;

%% Constants, functions and boundary conditions
N = 4;
xs = linspace(0, 1, N+2);
n = -3;

f = @(x) 0.*x;
a = @(x) 1 + x;

% f = @(x) exp(x);
% a = @(x) exp(x);

g = 1;

%% Analytical solution

% v = zeros(N+1, 1);
% for i = 1:length(xs)
%    v(i) = integral(@(y) (g + integral(@(z) f(z), y, 1)) / a(y), 0, xs(i), 'ArrayValued', true); 
% end
% 
% plot(xs, v)
% hold on

%% Calculations

integrand_A = @(x, i, j) a(x) * d_phi(i, x, xs) * d_phi(j, x, xs);
integrand_b = @(x, i) f(x) * phi(i, x, xs);

tic
A = calc_A(xs, integrand_A, N, n);
A(N+1, N+1) = quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N+1);
A(N+1, N) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N);
A(N, N+1) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N, N+1)

b = calc_b(xs, integrand_b, N, n, g);
b(N+1) = g;

c = [0; A\b];
toc

plot(xs, c)

%% Functions

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

function value = d_phi(i, x, xs)
    if i+2 > length(xs)
        i = i-1;
    end
    if x <= xs(i) || x >= xs(i+2)
        value = 0;
    elseif x >= xs(i) && x <= xs(i+1)
        value = 1/(xs(i+1) - xs(i));
    elseif x >= xs(i+1) && x <= xs(i+2)
        value = -1/(xs(i+2) - xs(i+1));
    end
end

function b = calc_b(xs, fun, N, n, g)
    b = zeros(N, 1);
    
    for i = 1:N
        b(i) = g*phi(i, 1, xs) + quadr_b(xs(i), xs(i+1), n, fun, i) + quadr_b(xs(i+1), xs(i+2), n, fun, i);
    end
end

function value = quadr_b(a, b, n, fun, i)
    if n == 1
        value = (b-a)*fun((b+a)/2, i);
    else    
        value = ((b-a)/2)*fun(((b-a)/2)*(-1/sqrt(3))+(b+a)/2, i) + ((b-a)/2)*fun(((b-a)/2)*(1/sqrt(3))+(b+a)/2, i);
    end
end

function A = calc_A(xs, fun, N, n)
    A = zeros(N, N);

    for i = 1:N
        for j = 1:N
            A(j, i) = quadr_A(xs(i), xs(i+1), n, fun, i, j) + quadr_A(xs(i+1), xs(i+2), n, fun, i, j);
        end
    end
end

function value = quadr_A(a, b ,n, fun, i, j)
    if n == 1
        value = (b-a)*fun((b+a)/2, i, j);
    else    
        value = ((b-a)/2)*fun(((b-a)/2)*(-1/sqrt(3))+(b+a)/2, i, j) + ((b-a)/2)*fun(((b-a)/2)*(1/sqrt(3))+(b+a)/2, i, j);
    end
end


