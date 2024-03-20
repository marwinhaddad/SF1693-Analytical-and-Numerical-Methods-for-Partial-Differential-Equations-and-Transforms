%% Initialvärden och konstanter

M = 1000;
N = 6;
K = 2*pi;

n = -N:N;

% --- diskretisering av området omega; Mxp är 1000 rader med varje punkt xp

Mrand = rand([M, 2]);


% --- skapa lista med alla n1 och n2 på samma sätt som Mxp
    %Mn 169x2 med varje kombination n1 n2
Mn = zeros((2*N+1)^2, 2);
pos = 1;

for m = 1:length(n)
    for l = 1:length(n)
        Mn(pos, :) = [n(m) n(l)];
        pos = pos + 1;
    end
end


% % --- skapa matrisen MX med fourierbas
MX = zeros(length(Mrand), length(Mn));
vx = [Mrand(:,2), 1-Mrand(:,1)];
fx = zeros(length(Mrand), 1);

for m = 1:length(Mrand)
    for l = 1:length(Mn)
        fx(m) = f(Mrand(m, :));
        MX(m, l) = exp(K*1i*(Mn(l, 1)*Mrand(m, 1) + Mn(l, 2)*Mrand(m, 2)));
    end
end

vhat1 = MX\vx(:,1);
vhat2 = MX\vx(:,2);
fhat = MX\fx;

size(vhat1)
