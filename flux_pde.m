
%% Problem Definition
% The following variables will define our problem:
%%
% * |g|: A specification function that is used by |initmesh|. For more
% information, please see the documentation page for |squarereg| and |pdegeom|.
% * |b|: A boundary file used by |assempde|. For more
% information on this file, please see the documentation pages for
% |squareb4| and |pdebound|.
% * |c|, |a|, |f|: The coefficients of the PDE.
g = @squareg;
c = -1;
a = 2;
f = '0.5*(y-1) + 2*u.^3'; % y0 = 1

bc = @(x) (1+tanh(x*5))/2;

%% Boundary Conditions
% Plot the geometry and display the edge labels for use in the boundary
% condition definition.
figure;
pdegplot(g, 'edgeLabels', 'on');
axis([-1.1 1.1 -1.1 1.1]);
title 'Geometry With Edge Labels Displayed'

% Create a pde entity for a PDE with a single dependent variable
pb = pde(1);
% Create a geometry entity
pg = pdeGeometryFromEdges(g);
innerBCFunc = @(thePde, loc, state) bc(loc.x);
b1 = pdeBoundaryConditions(pg.Edges(1), 'u', innerBCFunc); % north
b2 = pdeBoundaryConditions(pg.Edges(2), 'q', 0, 'g', 0); % west
b3 = pdeBoundaryConditions(pg.Edges(3), 'q', 0, 'g', 0); % south
b4 = pdeBoundaryConditions(pg.Edges(4), 'q', 0, 'g', 0); % east
pb.BoundaryConditions = [b1 b2 b3 b4];

%% Mesh for Fast Poisson Solver
%% Solve Using Both Fast Poisson Solver and Standard Solver
% We solve using both the fast Poisson solver that is implemented in
% |poisolv| and  the standard solver that is implemented in |assempde|.
n = 256;
tic;
u = pdenonlin(pb,p,e,t,c,a,f,'jacobian','lumped','report', 'on', ...
              'Tol', 4e-4);
tstandard = toc;
fprintf('%-5d|%15.5g\n',n,tstandard);

% calculate velocities
[ux,uy] = pdegrad(p,t,u);
U = -uy; V = ux;

%% plot solution
hf = figure;
insertAnnotation('flux_pde.m');
xvec = linspace(-1,1,50);
ax(4) = subplot(221);
plot(xvec, bc(xvec));
title(['f = ' f ' | BC'])
xlabel('x'); ylabel('\Psi(y = 1)');
ax(1) = subplot(2,2,3);
pdesurf(p,t,u);
colorbar('Location', 'SouthOutside');
caxis([-1 1]);
view(2)
xlabel('x');ylabel('y');title('\Psi');
liney(0.4);
figure(hf);
ax(2) = subplot(222);
pdesurf(p,t,U); view(2);
center_colorbar;
xlabel('x');ylabel('y');title('U');
ax(3) = subplot(224);
pdesurf(p,t,V); view(2);
center_colorbar;
xlabel('x');ylabel('y');title('V');
linkaxes(ax(1:3), 'y');
linkaxes(ax, 'x');