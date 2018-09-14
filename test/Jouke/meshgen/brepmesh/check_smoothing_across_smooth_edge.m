clc; clear; close all

% mesh
tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
idstyp = importdata('idstyp.dat');
v2h = importdata('mv2h.dat');
twin = importdata('mtwin.dat');

nt = size(tri,1);
nv = size(xyz,1);

% brep
verts = importdata('../../brep/verts.dat');
vxyz = verts(:,1:3);
edges = importdata('../../brep/edges.dat'); % face, orig, prev(2), next(2)
% faces = importdata('../../brep/faces.dat');
nverts = 0;
nedges = 0;
nfaces = 0;
for i = 1:nv
    switch idstyp(i,2)
        case 2
            nfaces = max(nfaces, idstyp(i,1));
        case 1
            nedges = max(nedges, idstyp(i,1));
        case 0
            nverts = max(nverts, idstyp(i,1));
    end
end




figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');


% l = [32668       32905];
% d = [-1.25166634E-03   2.77692801E-04   9.02278553E-05];
% x = [0.990000010      0.281401873      0.201836541];
% i = 287;

l = [31670       33163];
d = [1.05600478E-03  -1.29371893E-03  -4.20354772E-04];
x = [0.709069014      0.278593689      0.200924113];
i = 279;

% l = [];
% d = [];
% x = [];
% i = ;

plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'rs');
plot3(x(1), x(2), x(3), 'yo');
quiver3(xyz(l(1),1), xyz(l(1),2), xyz(l(1),3), d(1), d(2), d(3), 0, 'b', 'maxheadsize', 1);




aabb = zeros(1,6);
aabb(1:2:5) = min(xyz(l,:),[],1);
aabb(2:2:6) = max(xyz(l,:),[],1);
rng = aabb(2:2:6) - aabb(1:2:5);
% rng = max(rng, 0.1*max(rng));
rng(:) = 10*max(rng);
aabb(1:2:5) = aabb(1:2:5) - rng;
aabb(2:2:6) = aabb(2:2:6) + rng;


ih = 2*i-1;
next = edges(ih,5:6);
jh = 2*(next(1) - 1) + next(2);
k = [edges(ih,2), edges(jh,2)];

plot3(vxyz(k,1), vxyz(k,2), vxyz(k,3), 'g*--', 'linewidth', 1);

cl = colorcet( 'I2', 'N', 2 );
for j = 1:2
    ih = 2*(i-1) + j;
    f = edges(ih,1);
    k = find(idstyp(:,1) == f & idstyp(:,2) == 2);
    plot3(xyz(k,1), xyz(k,2), xyz(k,3), '.', 'color', cl(j,:), 'markersize', 12);
end


% axis image vis3d
axis(aabb);
daspect([1,1,1]);
axis vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end