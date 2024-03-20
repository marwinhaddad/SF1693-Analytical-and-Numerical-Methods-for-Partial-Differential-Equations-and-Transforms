clear all;
close all;

N = 200;
a = 1.0;
T =10;
c = 2;
dx = a / N;
dt = dx / 2.0;
M = fix(T/ dt);

figure(1)
axis([0 1 -2 2])
u = zeros(N-1,M+1);
p = zeros(N-1,M+1);
A = zeros(N-1,N-1);
X = zeros(N-1,1);
writeobj = VideoWriter('outw.avi');
writeobj.FrameRate = 10;
open(writeobj);
nframe = M;
mov(1:nframe)= struct('cdata', [], 'colormap' ,[]);
set(gca, 'nextplot', 'replacechildren');

for n = 1:N-1
    X(n, 1) = n*dx;
    A(n, n) = -2*N^2;
    if n == N-1
        break
    end
A(n+1, n) = N^2;
A(n, n+1) = N^2;
end

% g = @(x) sin(2 * pi * x); % Dirichlet
g = @(x) cos(( 1 / (N-1)) * x * pi); % Neumann
u(:, 1) = g(dx:dx:dx * (N-1));
u = symplectic(u, p, N, dt, A, c, M);

for m=1:M
    plot(X, u (:, m+1) , 'b' , 'Linewidth', 1)
    mov(m) = getframe(gcf);
end

writeVideo(writeobj, mov);
close(writeobj);
open(writeobj);
writeVideo(writeobj, mov);
close(writeobj);

function u = symplectic(u, p, N, dt, A, c ,M)
    for m = 1 :M
        for n = 1 :N-1
            p(n, m+1) = p(n, m) + dt * c^2 * A(n, :) * u(:, m); 
            u(n , m+1) = u(n, m) + dt * (p(n, m+1));
        end
    end
end