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





figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

% l = [2        4119];
% d = [0.00000000       2.39154627E-03   9.52202827E-04];
% i = 1;

l = [28        9967];
d = [0.00000000      -1.72558194E-03   4.41166078E-04];
i = 4;

% l = [];
% d = [];
% i = ;

plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'r.');
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

plot3(vxyz(k,1), vxyz(k,2), vxyz(k,3), 'g*-', 'linewidth', 1.5);


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