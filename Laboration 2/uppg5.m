clear;
clc;
close;

N = 999;
n = 2;

epsilon = 0.1;
xs = [epsilon:(1-epsilon)/N:1];
u = @(x) log(x)/log(epsilon);
g = @(x) 1/(epsilon - 1).*x + 1/(1 - epsilon);
d_g = 1/(epsilon - 1);

integrand_A = @(x, i, j) x * d_phi(i, x, xs) * d_phi(j, x, xs);
integrand_b = @(x, i) d_g * phi(i, x, xs);

A = calc_A(xs, integrand_A, N-1, n);
b = calc_b(xs, integrand_b, N-1, n);
w = [0; A\b; 0];
uh = w + g(xs)';


plot(xs, uh, 'b')
xline(0);
yline(0);
hold on
plot(-xs, uh, 'b')
axis([-1.5 1.5 -0.5 1.5]);

figure(2)
b = linspace(0, 2*pi, 100);
[R, A2] = ndgrid(xs,b);
[x2,y2]=size(R);
Z = ones(x2, y2).*uh;
[X,Y,Z] = pol2cart(A2,R,Z);

mesh(X, Y, Z)






