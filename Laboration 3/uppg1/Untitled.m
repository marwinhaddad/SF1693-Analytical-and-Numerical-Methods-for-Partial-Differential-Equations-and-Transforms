close all
clc

global N h X F f tol

N = 500000;
h = 0.0004;
a = 0.5;
tol = 1e-3;

X = zeros(N+1, 4);
X(1,:) = [1 - a, 0, 0, sqrt((1 + a)/(1 - a))];

F = @(x) [x(3), x(4), -x(1)/(x(1)^2 + x(2)^2)^(3/2), -x(2)/(x(1)^2 + x(2)^2)^(3/2)];
f = @(x) [-x(1)/(x(1)^2 + x(2)^2)^(3/2), -x(2)/(x(1)^2 + x(2)^2)^(3/2)];

X_exp = euler_exp();
X_imp = euler_imp();
X_simp = euler_simp();
X_mid = midpoint();

E_exp = get_energy(X_exp);
E_imp = get_energy(X_imp);
E_simp = get_energy(X_simp);
E_mid = get_energy(X_mid);

figure(1)
plot(X_exp(:, 1), X_exp(:, 2))
axis equal
hold on
plot(X_imp(:, 1), X_imp(:, 2));
plot(X_mid(:, 1), X_mid(:, 2));
plot(X_simp(:, 1), X_simp(:, 2));

figure(2)
plot(1:N+1, E_exp)
hold on 
plot(1:N+1, E_imp)
plot(1:N+1, E_mid)
plot(1:N+1, E_simp)


