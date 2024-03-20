% Vagekvationen u``(t)=Au med Dirichlet randvillkor u(t,0)=u(t,1)=0

clear all
close all

N=200; % antal interval
a=1.0; % parameter i begynnelsedata
T=2; % sluttid

dx=a/(N+1); % steglangd i rummet
dt=dx/1; % tidssteg, tank på stabilitetsvillkoren
M=fix(T/dt);% antal tidsteg

% allokering av minne
u=zeros(N,M); % u(n,m) utbojning vid tid m*dt i position n*dx
p=zeros(N,M); % p=u’
A=zeros(N,N); % Au är differensapproximation av d^2 u/dx^2
X=zeros(N,1); % X(n) blir n*dx

vv = VideoWriter('outw.avi');
open(vv);

BV = @(x) sin(5*pi*x);
BVt = @(x) 0;

for n=1:N % slinga med rumssteg för att bilda begynnelsedata och A
    X(n) = n * dx;
    u(n, 1) = BV(X(n));
    p(n, 1) = BVt(X(n));
    
    A(n, n) = -2*N^2;
    if n == N
        break
    end
    A(n, n+1) = N^2;
    A(n+1, n) = N^2;
end

% randvillkor
u(1, :) = 0;
u(end, :) = 0;

% p(1, :) = 0;
% p(end, :) = 0;   

for m=1:M % tidstegning med symplektiska Euler
    u(:, m+1) = u(:, m) + dt .* p(:, m);
    p(:, m+1) = p(:, m) + dt .* (A * u(:, m+1));
    
    figure(1)
    plot(X,u(:,m), 'b', 'Linewidth', 1)
    axis([0 1 -1.5 1.5])
    frame = getframe(gcf);
    writeVideo(vv, frame);
end
close(vv);