clc; clear; close all

% mesh
tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
idstyp = importdata('idstyp.dat');
v2h = importdata('mv2h.dat');
twin = importdata('mtwin.dat');

xyz = importdata('xyz_smooth.dat');
idstyp = importdata('idstyp_smooth.dat');

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


% l = 1936;
% d = [6.6264624860738243E-004   1.5021649326273477E-005   7.7580044752927478E-005];
% x = [0.99270217539309979      -0.10495942556589022       0.25922652835553422];
% y = [0.99203970247832007      -0.10497688377466710       0.25914798402273181];
% q = [6.6264624860734589E-004   1.5021649326277253E-005   7.7580044752945110E-005];
% i = 1;
% h = [36           2];
% z = [
%     1.0006582472848098      -0.10436811035834220       0.25981874645354314
%     1.0006579779595062      -0.10436670244456078       0.25982126981844800
%     ];


l = 2045;
ie = 33;
2*(ie - 1) + [1,2]
faces = edges(2*(ie - 1) + [1,2],1)


return






l = 1989;
d = [-1.1444316968277186E-005  -2.3634044636028682E-004  -1.6344178466770040E-005];
x = [0.53352479534075647       -1.5236327155883300E-002  0.23325439083150182];
y = [0.53353632959835395       -1.4999999999999999E-002  0.23327086345849901];
q = [-1.1444316968312052E-005  -2.3634044636028709E-004  -1.6344178466770609E-005];
i = 1;
h = [17           2];
z = [
    0.53352488528138575       -1.5236331544999633E-002  0.23325439955803226     
    0.53352488528138575       -1.5236340446360287E-002  0.23325451928003224
    ];

% l = [];
% d = [];
% x = [];
% i = ;

plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'ro');
quiver3(xyz(l,1), xyz(l,2), xyz(l,3), d(1), d(2), d(3), 0, 'b', 'maxheadsize', 1);
plot3(x(1), x(2), x(3), 'y*');
plot3(y(1), y(2), y(3), 'go');
quiver3(y(1), y(2), y(3), q(1), q(2), q(3), 0, 'g', 'maxheadsize', 1);
plot3(z(:,1), z(:,2), z(:,3), 'cs', 'markersize', 10 );


% ih = (2*h(1)-1) + h(2);
% next = edges(ih,5:6);
% jh = 2*(next(1) - 1) + next(2);
% k = [edges(ih,2), edges(jh,2)];
% 
% plot3(vxyz(k,1), vxyz(k,2), vxyz(k,3), 'ms--', 'linewidth', 1);


ih = 2*(h(1) - 1) + h(2);
next = edges(ih,5:6);
kh = 2*(next(1) - 1) + next(2);
lv = [edges(ih,2), edges(kh,2)];
plot3( vxyz(lv,1), vxyz(lv,2), vxyz(lv,3), 'ms--', 'linewidth', 1);




aabb = zeros(1,6);
aabb(1:2:5) = xyz(l,:);
aabb(2:2:6) = xyz(l,:);
rng = 0.1*ones(1,3);
aabb(1:2:5) = aabb(1:2:5) - rng;
aabb(2:2:6) = aabb(2:2:6) + rng;


% ih = 2*i-1;
% next = edges(ih,5:6);
% jh = 2*(next(1) - 1) + next(2);
% k = [edges(ih,2), edges(jh,2)];
% 
% plot3(vxyz(k,1), vxyz(k,2), vxyz(k,3), 'g*--', 'linewidth', 1);

% cl = colorcet( 'I2', 'N', 2 );
% for j = 1:2
%     ih = 2*(i-1) + j;
%     f = edges(ih,1);
%     k = find(idstyp(:,1) == f & idstyp(:,2) == 2);
%     plot3(xyz(k,1), xyz(k,2), xyz(k,3), '.', 'color', cl(j,:), 'markersize', 12);
% end


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