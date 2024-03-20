clc
clear
close all;

N = 2000;

L = pi;
dx = L / (N+1);
dt = 0.01;

x = 0:dx:L;
t = 0:dt:L;

n = 10000;
nnot = n/200;

% --- skapa matrisen
M = zeros(length(x), length(t));
Mn = M;

% --- räkna ut begynnelsevärde
Mn(:, 1) = (pi - x).*x;

for nn = 1:n
    fun = @(x) (pi - x).*x.*sin(nn.*x); 
    an = 2/pi * integral(fun, 0, L);
    for tt = 2:length(t) - 1
        for xx = 2:length(x) - 1
            Mn(xx, tt) = Mn(xx, tt) + an .*sin(nn.*x(xx)).*exp(-nn^2.*t(tt));
        end
    end
end

medel_kvadrat_fel = zeros(nnot, 1);
M(:, 1) = (pi - x).*x;
for nn = 1:nnot
    fun = @(x) (pi - x).*x.*sin(nn.*x); 
    an = 2/pi * integral(fun, 0, L);
    for tt = 2:length(t) - 1
        for xx = 2:length(x) - 1
            M(xx, tt) = M(xx, tt) + an .*sin(nn.*x(xx)).*exp(-nn^2.*t(tt));
        end
    end
    medel_kvadrat_matris = abs(Mn-M).^2;
    medel_kvadrat_fel(nn, 1) = sum(medel_kvadrat_matris(:))/numel(Mn);
end

plot(medel_kvadrat_fel, 'LineWidth',2.0)
