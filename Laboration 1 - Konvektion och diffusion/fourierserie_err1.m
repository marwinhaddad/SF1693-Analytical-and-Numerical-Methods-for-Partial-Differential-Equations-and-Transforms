clc
clear
close all;

N = 400;

L = [0, pi]; % Intervall i x-led
dx = (L(2) - L(1))/ (N+1);
x = L(1):dx:L(2);

I = [0, pi]; % Intervall i t-led
dt = (I(2) - I(1))/ (N+1);
t = I(1):dt:I(2);

n = 100;
% 
% % --- skapa matrisen
% 
% M = zeros(length(x), length(t));
% 
% % --- räkna ut begynnelsevärde när t = 0
% 
% M(:, 1) = (pi - x) .* x;
% 
% % --- räkna varje a_n
% a_n = zeros(n, 1);
% 
% for nn = 1:n
%     fun = @(s) (pi - s) .* s .* sin(nn .* s);
%     a_n(nn , 1) = (2/pi) * integral(fun, 0, pi);
% end
% 
% for nn = 1:n
%     for jj = 2:length(t) - 1
%         for ii = 2:length(x) - 1
%            M(ii, jj) = M(ii, jj)+a_n(nn)*sin(nn*x(ii))*exp(-nn^2*t(jj)); 
%         end
%     end
% end



