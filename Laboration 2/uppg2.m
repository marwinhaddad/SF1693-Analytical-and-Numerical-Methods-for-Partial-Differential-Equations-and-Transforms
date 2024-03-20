clear all;
clc;
close all;

%% Constants and functions
N = 9;  % Inner discretization points
xs = linspace(0, 1, N+2);
n = 1;

f = @(x) 0.*x;
a = @(x) 1 + x;
g = 1;

%% Analytical solution
for i = 1:length(xs)-1
   v(i) = integral(@(y) (g + integral(@(z) f(z), y, 1)) / a(y), 0, xs(i+1), 'ArrayValued', true); 
end

%% Calculations

integrand_A = @(x, i, j) a(x) * d_phi(i, x, xs) * d_phi(j, x, xs);
integrand_b = @(x, i) f(x) * phi(i, x, xs);

A = calc_A(xs, integrand_A, N, n);
A(N+1, N+1) = quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N+1);
A(N+1, N) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N);
A(N, N+1) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N, N+1);

b = calc_b(xs, integrand_b, N, n);
b(N+1) = g;

c = [A\b];

err = c - v;
E_norm = norm(err, inf)

























