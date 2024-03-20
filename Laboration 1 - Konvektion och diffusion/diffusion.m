clc
clear
close all;
format long

%% Initialvärden och konstanter

N = 100;

L = [0, pi]; % Intervall i x-led
dx = (L(2) - L(1))/ (N+1);
x = L(1):dx:L(2);

I = [0, pi]; % Intervall i t-led
dt = (dx^2)/2;
t = I(1):dt:I(2);
disp(dt/dx^2)
h = dt/dx^2;
k = 1 - 2*h;

%% Numerisk lösning FD: g(x) = x(pi - x)
M = zeros(length(x), length(t));
M(:, 1) = (-x + pi).*x;

for j = 1:length(t) - 1
    for i = 2:length(x) - 1
       M(i, j+1) = h*M(i-1, j) - k*M(i, j) + h*M(i+1, j);
    end 
end
% --- Plot
figure
[tt, xx] = meshgrid(t, x);
mesh(xx, tt, M)
xlabel('Längdenhet')
ylabel('Tidsenhet')
zlabel('Temperatur')
axis([0 pi 0 pi])


%% Nummerisk lösning FD: g(x) = { x ; x = [0, pi/2), 0 ; x = [pi/2, pi]}
M = zeros(length(x), length(t));

for i = 1:length(x)
    if x(i) < pi/2
        M(i, 1) = x(i);
    else
        M(i, 1) = 0;
    end
end

for j = 1:length(t) - 1
    for i = 2:length(x) - 1
       M(i, j+1) = h*M(i-1, j) - k*M(i, j) + h*M(i+1, j);
    end 
end
% --- Plot
figure
[tt, xx] = meshgrid(t, x);
mesh(xx, tt, M)
xlabel('Längdenhet')
ylabel('Tidsenhet')
zlabel('Temperatur')
axis([0 pi 0 pi])


%% Error Finita Differenser

% Err = zeros(length(x), length(t));
% 
% fun = @(x) (pi - x).*x.*sin(x);
% a1 = (2/pi) * integral(fun, 0, pi);
% for j = 1:length(t) - 1
%     for i = 2:length(x) - 1
%         Err(i, j+1) = sin(x(i))*exp((-1^2)*j);
%     end
% end
% E = abs(Err - M);
% figure()
% contour(E, 200, 'linecolor', 'non')
% colormap(jet(256))
% colorbar
