clear;
clc;
close;

N = 99;
n = 2;

epsilon = 0.0;
xs = [epsilon:(1-epsilon)/N:1];
u = @(x) log(x)/log(epsilon);
g = @(x) 1/(epsilon - 1).*x + 1/(1 - epsilon);
d_g = 1/(epsilon - 1);

integrand_A = @(x, i, j) x * d_phi(i, x, xs) * d_phi(j, x, xs);
integrand_b = @(x, i) d_g * phi(i, x, xs);

A = calc_A(xs, integrand_A, N-1, n);
b = calc_b(xs, integrand_b, N-1, n);
w = [0; A\b; 0];
uh = w + g(xs)'

figure(1)
plot(xs, uh, 'b')
xline(0);
yline(0);
hold on
plot(-xs, uh, 'b')
axis([-1.5 1.5 -0.5 1.5]);

%% nedan kod tagen fran http://math.loyola.edu/~loberbro/matlab/html/PlotChangeCoordinates.html
figure(2)
b = linspace(0, 2*pi, 100);
[R, A2] = ndgrid(xs,b);
[x2,y2]=size(R);
Z = ones(x2, y2).*uh;
[X,Y,Z] = pol2cart(A2,R,Z);
colormap(jet(256))
mesh(X, Y, Z)

