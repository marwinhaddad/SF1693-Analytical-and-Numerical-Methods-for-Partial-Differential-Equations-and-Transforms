close all
clc

N = 500000;
h = 0.0004;
a = 0.5;

Tint = [0:h:N*h];

BV = [1 - a; 0; 0; sqrt((1 + a)/(1 - a))];
F = @(t, X) [X(3); X(4); -X(1)/(X(1)^2 + X(2)^2)^(3/2); -X(2)/(X(1)^2 + X(2)^2)^(3/2)];

[t, x] = ode45(F, Tint , BV);

scatter(x(1, 1), x(1,2),'g', 'filled')
hold on
scatter(x(end, 1), x(end, 2), 'r', 'filled')
scatter(0,0, 'black', 'filled')
plot(x(:, 1), x(:, 2),'black', 'LineWidth', 0.001);
axis equal
legend('Start', 'Stop', 'Origo', 'Orbit')
