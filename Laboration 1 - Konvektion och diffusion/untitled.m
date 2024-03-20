N = 10;

Lx = 12;
Ly = 5;

dx = Lx/N;
dy = Ly/N;

x = 0:dx:Lx;
y = 0:dy:Ly;

Text = 25;

T = zeros(length(x), length(y));

% randvillkor
T(end, :) = Text;


% alla andra randvillkor är noll neumann


for i = 2:length(x) - 1
    for j = 2:length(y) - 1
        T(i, j) = 1/4 * (T(i-1, j) + T(i+1, j) + T(i, j-1) + T(i, j+1));
    end
end

figure
[xx, yy] = meshgrid(x, y);
mesh(xx, yy, T)
