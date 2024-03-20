clear;
clc;
close;

I_N = [3, 99]; 

for N = I_N(1):I_N(2)
    xs = linspace(0, 1, N+2);
    h = xs(2) - xs(1);
    n = 1;

    f = @(x) exp(x);
    a = @(x) exp(x);
    g = 1;

    integrand_A = @(x, i, j) a(x) * d_phi(i, x, xs) * d_phi(j, x, xs);
    integrand_b = @(x, i) f(x) * phi(i, x, xs);
    
    A = calc_A(xs, integrand_A, N, n);
    A(N+1, N+1) = quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N+1);
    A(N+1, N) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N+1, N);
    A(N, N+1) = - quadr_A(xs(end-1), xs(end) ,n, integrand_A, N, N+1);

    b = calc_b(xs, integrand_b, N, n);
    b(N+1) = g;
    tic, x = A\b; t = toc;
    U_h = [0; x];
    
    U = zeros(length(xs), 1); dU = U(1:end-1); dU_h = U(1:end-1);

    for i = 1:length(xs)
       U(i) = integral(@(y) (g + integral(@(z) f(z), y, 1)) / a(y), 0, xs(i), 'ArrayValued', true);
       if i == length(xs)
           i = i - 1;
       end
       dU(i) = (g + integral(@(x) f(x), xs(i), 1))/ a(xs(i));
       dU_h(i) = (U_h(i+1) - U_h(i))/h;
    end

    E_norm = 0;
    for i = 1:length(xs)-1
        E_norm = E_norm + integral(@(x) a(x) * (dU(i) - dU_h(i))^2, xs(i), xs(i+1));    
    end
    t = toc;
    E_vec(N-(I_N(1)-1)) = E_norm;
    t_vec(N-(I_N(1)-1)) = t;
end
figure(1)
plot(I_N(1):I_N(2), E_vec)

figure(2)
plot(I_N(1):I_N(2), t_vec)
