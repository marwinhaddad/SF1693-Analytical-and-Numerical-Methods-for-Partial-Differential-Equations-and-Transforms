format long
%a)
L = 3;
k = 2;
T0 = 290; Tn1 = 400; 
q0 = 3000; q1 = 200;
n = 5;
h = (L/(n+1));
x0 = 0; xn = L;
x = [x0+h:h:xn-h]';
Q = @(x) q0.*exp(-q1.*(x-0.7.*L).^2)+200;

A = (diag(2*k*ones(n,1)) + diag(-1*k*ones(n-1,1),1) + diag(-1*k*ones(n-1,1),-1))*(1/h^2)

b = zeros(n, 1);
b(1) = Q(x(1))+T0*(k/h^2);
b(2:n-1) = Q(x(2:n-1));
b(n) = Q(x(end))+Tn1*(k/h^2);


%b)
T = A\b;
T = [T0; T; Tn1]
plot(T)

%c)
nVec = zeros(4, 1);
Tmax = zeros(4, 1); Tmin = zeros(4, 1); Tmed = zeros(4, 1);
n = 40;
for jj = 1:4
    h = (L/(n+1));
    x0 = 0; xn = L;
    x = [x0+h:h:xn-h]';
    
    A = (diag(2*ones(n,1)) + diag(-1*ones(n-1,1),1) + diag(-1*ones(n-1,1),-1))*(k/h^2);
    
    b = zeros(n, 1);
    b(1) = Q(x(1))+T0*(k/h^2);
    b(2:n-1) = Q(x(2:n-1));
    b(n) = Q(x(end))+Tn1*(k/h^2);
    
    T = A\b;
    T = [T0; T; Tn1];
    plot([0;x;3],T)
    hold on
    Tmax(jj) = max(T);
    Tmin(jj) = min(T);
    Tmed(jj) = mean(T);
    nVec(jj) = n;
    n = n*2;
end
legend('n = 40', 'n = 80', 'n = 160', 'n = 320')
tab = table(nVec, Tmax, Tmin, Tmed);
tab.Properties.VariableNames = {'n-Värden', 'Max', 'Min', 'Medel'}
cond(A,inf)
% Vi ser av konditionstalet för matris A att talet är väldigt stort, detta
% innebär att oavsett hur precisa våra parametrar q1 och q0 är så kommer vi
% få ett fel i vårt svar

