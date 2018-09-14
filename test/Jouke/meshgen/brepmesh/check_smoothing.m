clc; clear; close all

% mesh
tri = importdata('tri.dat');
xyz = importdata('xyz.dat');
idstyp = importdata('idstyp.dat');

xyzs = importdata('xyz_smooth.dat');

% xyz = importdata('xyz_smooth.dat');
% idstyp = importdata('idstyp_smooth.dat');

% vertex displacement
dxyz = importdata('debug_dxyz.dat');




figure('units','normalized','position',[.15 .15 .7 .7]);
hold on

trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ...
    'facecolor', 'w', 'specularstrength', 0, 'edgecolor', 'k');

% quiver3(xyz(:,1), xyz(:,2), xyz(:,3), ...
%     dxyz(:,1), dxyz(:,2), dxyz(:,3), ...
%     0, 'color', 'r', 'maxheadsize', 1);

% l = 4243;
% l = 19408;
l = 4274;
plot3( xyz(l,1), xyz(l,2), xyz(l,3), 'g*');%'r*' );
quiver3(xyz(l,1), xyz(l,2), xyz(l,3), ...
    dxyz(l,1), dxyz(l,2), dxyz(l,3), ...
    0, 'color', 'r', 'maxheadsize', 1);
% k = 4261;
% k = [19409, 4261];
k = [1936, 4261];
% plot3( xyz(k,1), xyz(k,2), xyz(k,3), 'g*' ); 

axis image vis3d
view(3)
camproj('persp');
if 1
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end


idstyp(l,:)
idstyp(k,:)


norm( xyz(l,:) - xyzs(l,:) )