clear;
clc;
close;

N = 99;
epsilon = .1;
xs = [epsilon:(1 - epsilon)/(N):1];
h = xs(2) - xs(1);
sigma = 0.001;
TOL = 0.01;
diff = 1;

integrand_A = @(x, i ,j) x * d_phi(i, x, xs) * d_phi(j, x, xs);
integrand_b = @(x, i) 1/(epsilon - 1) * phi(i, x, xs);

A = calc_A(xs, integrand_A, N-1, 2);
u = [1; ones(N-1, 1).*0.5; 0];

A(N, N) = 1/(2*h^2) * (xs(end)^2 - xs(end-1)^2);
A(N, N-1) = -A(N, N);
A(N-1, N) = A(N, N-1);

A = [zeros(N, 1), A];
A = [zeros(1, N+1); A];

A(1, 1) = 1/(2*h^2) * (xs(2)^2 - xs(1)^2);
A(1, 2) = -A(1, 1);
A(2, 1) = A(1, 2);
for n = 1:N+1
    for k = 2:N
        u(k, n+1) = u(k, n) - sigma * dF(A, u(:,n), k);
    end
end

u(1, end) = 1
plot(xs, u(:,end))

function value = dF(A, u, k)
    value = 0;
    for i = 1:length(u)
        value = value + A(i, k)*u(i);
    end
end



