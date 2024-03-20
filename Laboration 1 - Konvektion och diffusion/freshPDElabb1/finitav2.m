clc
clear
close all;

N = 100;

L = pi;
dx = L / (N+1);
dt = dx^2/2;

h = dt/dx^2;
k = 1 - 2*h;

x = 0:dx:L;
t = 0:dt:L;

% --- skapa matrisen
M = zeros(length(x), length(t));

% --- räkna ut begynnelsevärde
m = 1;
while x(m) < pi/2
    M(m, 1) = x(m);
    m = m + 1;
end

% --- explicit euler med finita differenser
for tt = 1:length(t) - 1
    for xx = 2:length(x) - 1
        M(xx, tt+1) = h * M(xx-1, tt) - k * M(xx, tt) + h * M(xx+1, tt);
    end
end

[T, X] = meshgrid(t, x);
mesh(X, T, M)
colormap(jet(256))
xlabel('Längdenhet')
ylabel('Tidsenehet')
zlabel('Temperatur')
axis([0 pi 0 pi])