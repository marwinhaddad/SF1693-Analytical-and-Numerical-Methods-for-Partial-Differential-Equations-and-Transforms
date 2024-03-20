clear all;
clc;
close all;

%% Constants and functions
N = 4;
xs = linspace(0, 1, N+2);

for n = 2

    f = @(x) 0;
    a = @(x) 5 - 4*x;
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

    A = calc_A(xs, integrand_A, N, n);
    A(N+1, N+1) = quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N+1);
    A(N+1, N) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N);
    A(N, N+1) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N, N+1);
    
    disp(A)

    b = calc_b(xs, integrand_b, N, n);
    b(N+1) = g;

    c = [0; A\b];

    plot(xs, c)
    hold on
end
