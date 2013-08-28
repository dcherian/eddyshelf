%% create

center = [20,30,40];
radii  = [40,20,10] .* [1000 1000 1]; % (x,y) in km, z in m
N = (10); % total number of floats = N^3

[fx,fy,fz] = ellipsoid(center(1),center(2),center(3),radii(1),radii(2),radii(3),N);

surf(fx/1000,fy/1000,fz);

%% read data & find convex hull

P = gallery('uniformdata',[25,3],1);
DT = delaunayTriangulation(P);
[K,vol] = convexHull(DT);

% plot hull
trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),...
       'FaceColor','cyan')