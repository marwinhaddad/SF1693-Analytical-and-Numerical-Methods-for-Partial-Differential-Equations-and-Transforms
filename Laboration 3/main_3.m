% Vågekvationen u’’(t)=Au med Dirichlet randvillkor u(t,0)=u(t,1)=0
clear all
close all

N = 10; % antal interval
a = 1.0; % parameter i begynnelsedata
T = 1.0; % sluttid
dx = a/(N+1); % steglängd i rummet
dt =dx/1.1; % tidssteg, tänk på stabilitetsvillkoren

c = 1;

M = fix(T/dt); % antal tidsteg

figure(1)
axis([0 1 -1 1])

% allokering av minne
u = zeros(N, M); % u(n,m) utböjning vid tid m*dt i position n*dx
p = zeros(N, M); % p=u 
A = zeros(N, N); % Au är differensapproximation av d^2 u/dx^2
X = zeros(N, 1); % X(n) blir n*dx
v = zeros(N, M);
vv = VideoWriter('outw.avi');
open(vv);

for n = 1:N % slinga med rumssteg för att bilda begynnelsedata och A

    X(n) = n * dx;
    % Dirichlet
    u(n, 1) = sin((n*pi)/a * X(n));

    % Neumann
    % u(n, 1) = cos(pi * X(n));
    
    if n == N
        break
    end
    
    A(n, n) = -2 ;
    A(n+1, n) = 1;
    A(n, n+1) = 1;
end

% bilda A på randen
% A(1, :) = 0;
% A(end, :) = 0;

for m = 1:M % tidstegning med symplektiska Euler
    for n = 1:N
        p(n, m+1) = p(n, m) + X(n) .* (c^2 .* N^2 .* (A(n, :) * u(:, m)));
        u(n, m+1) = u(n ,m) + X(n) .* p(n, m+1);
    end

    plot(X, u(:, m), 'b', 'Linewidth', 1)
    % axis([])
    frame = getframe(gcf);
    writeVideo(vv, frame);
    
end

close(vv);