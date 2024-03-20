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
M(:, 1) = (pi - x).*x;
W = M;

% --- explicit euler med finita differenser

for antalhalveringar = 1:10
    for tt = 1:length(t) - 1
        for xx = 2:length(x) - 1
            M(length(x)/2, tt+1) = dt/dx^2 * M(length(x)/2-1, tt) - 1 - 2 * dt/dx^2 * M(length(x)/2, tt) + dt/dx^2 * M(length(x)/2+1, tt);
        end
    end
    dt = dt / 2;
    for tt = 1:length(t) - 1
        for xx = 2:length(x) - 1
            W(length(x)/2, tt+1) = dt/dx^2 * W(length(x)/2-1, tt) - 1 - 2 * dt/dx^2 * W(length(x)/2, tt) + dt/dx^2 * W(length(x)/2+1, tt);
        end
    end
    halveringar(antalhalveringar) = abs(max(M(length(x)/2, :) - W(length(x)/2, :)));
end

for l = 1:length(halveringar)-1
    P(l) = halveringar(l)/halveringar(l+1);
end
halveringar  = halveringar(1:length(P));
tab = table(halveringar', P', log10(P')./log10(2));
tab.Properties.VariableNames = {'Fel', 'Kvoter', 'Noggrannhetsordningar'};
disp(tab)
