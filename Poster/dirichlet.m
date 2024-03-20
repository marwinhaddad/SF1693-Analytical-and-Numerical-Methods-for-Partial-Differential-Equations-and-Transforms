clc
close all

k = 400; % Specifik värmeledningsförmåga [W/(m * k)]
c = 386; % Specifik värmekapacitet [J/(kg * K)]

rho = 8960; % Densitet [kg/m^3]
y = 1; % Höjd [m]
x = 1; % Bredd [m]
z = 0.01; % Tjocklek [m]

boltz = 5.670374419e-8; % Stefan-Boltzmanns konstant [W/(m^2 * K^4]
h = 1; % Värmeldeningscoefficient [W/ (m^2 * K)]
Tsur = 300; % Omgivningens temperatur antas vara 293.15 K (20 C)
Tplate = 1000; % [K]
em = 0.5; % Förmåga att avge infraröd energi

PDEcount = 1;
model = createpde(PDEcount);

square = [3 4 0 x x 0 0 0 y y]';
g = decsg(square, 'S1', ('S1')');

geometryFromEdges(model, g);

figure;
pdegplot(model, 'EdgeLabels', 'on');
axis([-0.1 1.1 -0.1 1.1]);
title('Geometry With Edge Labels Displayed');


kz = z * k; % Thermic inductance

a = @(~, state) 2 * h + 2 * em * boltz * state.u.^3;
f = 2 * h * Tsur + 2 * em * boltz * Tsur^4;
d = z * rho * c;

specifyCoefficients(model, 'm', 0, 'd', 0, 'c', kz, 'a', a, 'f', f);

applyBoundaryCondition(model, 'dirichlet', 'Edge', 1, 'u', Tplate);
setInitialConditions(model, 0);

hmax = 0.1; % storlek på trianglar
M = generateMesh(model, 'Hmax', hmax);

figure;
pdeplot(model);
axis([-0.1 1.1 -0.1 1.1]);
title 'Plate With Triangular Elements Mesh'
xlabel 'X [m]'
ylabel 'Y [m]'

R = solvepde(model);
u = R.NodalSolution;

figure;
pdeplot(model, 'XYData', u, 'Contour', 'on', 'ColorMap', 'jet');
title 'Temperature In The Plate, Steady State Solution'
xlabel 'X [m]'
ylabel 'Y [m]'
axis([-0.1 1.1 -0.1 1.1]);


p = M.Nodes;

plotAlongY(p, u, 0);
title 'Temperature Along the Hight of Plate'
xlabel 'Y [m]'
ylabel 'T [K]'


specifyCoefficients(model, 'm', 0, 'd', d, 'c', kz, 'a', a, 'f', f);
tEnd = 5000;
tVec = 0:50:tEnd;
numNodes = size(p, 2);

u0(1:numNodes) = 300;

setInitialConditions(model, Tplate, 'edge', 1);

model.SolverOptions.RelativeTolerance = 1e-3;
model.SolverOptions.AbsoluteTolerance = 1e-4;

R = solvepde(model, tVec);
u = R.NodalSolution;

figure;
plot(tVec, u(3, :));
grid on
title 'Temperature Along the Hight of Plate as Function of Time'
xlabel 'Time [s]'
ylabel 'T [K]'

figure;
pdeplot(model, 'XYData', u(:, end), 'Contour', 'on', 'ColorMap', 'jet');
title(sprintf('Temperature of Plate, Transient Solution (%d s)', tVec(1,end)));
xlabel 'X [m]'
ylabel 'Y [m]'









