

N = 10;
a = @(x) 1 + x;
f = @(x) 0.*x;
h = 0.1;
g = 1;

v = zeros(N,1);
x = [0:h:1]';
for i = 1:length(x)
    v(i) = integral(@(y) (g + integral(@(z) f(z), y, 1)) / a(y), 0, x(i), 'ArrayValued',true);
end
plot(x, v)