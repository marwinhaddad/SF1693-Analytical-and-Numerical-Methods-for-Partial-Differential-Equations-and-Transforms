clc
clear
close all;

N = 1000;

L = pi;
dx = L / (N+1);
dt = 0.01;

x = 0:dx:L;
t = 0:dt:L;

n = 100;

% --- skapa matrisen
M = zeros(length(x), length(t));

% --- räkna ut begynnelsevärde
M(:, 1) = (pi - x).*x;

for nn = 1:n
    fun = @(x) (pi - x).*x.*sin(nn.*x); 
    an = 2/pi * integral(fun, 0, L);
    for tt = 2:length(t) - 1
        for xx = 2:length(x) - 1
            M(xx, tt) = M(xx, tt) + an .*sin(nn.*x(xx)).*exp(-nn^2.*t(tt));
        end
    end
end

[T, X] = meshgrid(t, x);
mesh(X, T, M)
colormap(jet(256))
xlabel('Längdenhet')
ylabel('Tidsenehet')
zlabel('Temperatur')
axis([0 pi 0 pi])