close all
clear all
clc

global N h X F f tol

N = 500000;
h = 0.0004;
a = 0.5;
tol = 1e-3;

X = zeros(N+1, 4);

% X(.,:) = [q1, q2, p1, p2]
X(1,:) = [1 - a, 0, 0, sqrt((1 + a)/(1 - a))];

F = @(x) [x(3), x(4), -x(1)/(x(1)^2 + x(2)^2)^(3/2), -x(2)/(x(1)^2 + x(2)^2)^(3/2)];
f = @(x) [-x(1)/(x(1)^2 + x(2)^2)^(3/2), -x(2)/(x(1)^2 + x(2)^2)^(3/2)];

X_exp = euler_exp();
% X_imp = euler_imp();
% X_symp = euler_symp();
% X_mid = midpoint();

E_exp = get_energy(X_exp);
% E_imp = get_energy(X_imp);
% E_symp = get_energy(X_symp);
% E_mid = get_energy(X_mid);

figure(1)
axis equal
hold on
scatter(X_exp(1, 1), X_exp(1, 2),'g','filled');
scatter(X_exp(end, 1), X_exp(end, 2), 'r', 'filled')
scatter(0,0, 'black', 'filled')
plot(X_exp(:, 1), X_exp(:, 2), 'black', 'LineWidth', 0.0001);
legend('Start', 'Stop', 'Origo', 'Orbit')
hold off
% 
% figure(2)
% axis equal
% hold on
% scatter(X_imp(1, 1), X_imp(1, 2),'g','filled');
% scatter(X_imp(end, 1), X_imp(end, 2), 'r', 'filled')
% scatter(0,0, 'black', 'filled')
% plot(X_imp(:, 1), X_imp(:, 2), 'black', 'LineWidth', 0.0001);
% legend('Start', 'Stop', 'Origo', 'Orbit')
% hold off
% 
% figure(3)
% axis equal
% hold on
% scatter(X_mid(1, 1), X_mid(1, 2),'g','filled');
% scatter(X_mid(end, 1), X_mid(end, 2), 'r', 'filled')
% scatter(0,0, 'black', 'filled')
% plot(X_mid(:, 1), X_mid(:, 2), 'black', 'LineWidth', 0.0001); 
% legend('Start', 'Stop', 'Origo', 'Orbit')
% hold off
% 
% figure(4)
% axis equal
% hold on
% scatter(X_symp(1, 1), X_symp(1, 2),'g','filled');
% scatter(X_symp(end, 1), X_symp(end, 2), 'r', 'filled')
% scatter(0,0, 'black', 'filled')
% plot(X_symp(:, 1), X_symp(:, 2), 'black', 'LineWidth', 0.0001);
% legend('Start', 'Stop', 'Origo', 'Orbit')
% hold off

figure(5)
plot(1:N+1, E_exp)
% 
% 
% figure(6)
% plot(1:N+1, E_imp)
% 
% figure(7)
% plot(1:N+1, E_mid)
% 
% figure(8)
% plot(1:N+1, E_symp)
% 
% figure(9)
% plot(1:N+1, E_exp)
% hold on
% plot(1:N+1, E_imp)
% plot(1:N+1, E_mid)
% plot(1:N+1, E_symp)
% legend('Explicit', 'Implicit', 'Mid point', 'Symplectic', 'Location','southwest')


