clc; clear; close all

tri = importdata('../connect_01.dat');
xyz = importdata('../pos_000.dat');

fprintf('%d verts, %d triangles\n',size(xyz,1), size(tri,1));

l = 84;%[1666,2773];
d = [
    0.0000000000000000       -4.4064128847254769E-005  -6.4556588239085653E-004
    ];
% -3.7570408007703004E-018  -5.1616942450078231E-004   5.1616942450104176E-004];

p = [
    -2.9999999999999999E-002 -0.14894395385626549        1.7767909546782489E-002
  -2.9999999999999999E-002 -0.14979673993204731        7.8061966238738151E-003
%     0.25505503005977892        2.4748737341529083E-002   2.4748737341529173E-002
%   0.25505503005977892       -2.4748737341529083E-002   2.4748737341529173E-002
];

figure;
hold on
trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'w', 'specularstrength', 0);

% plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'r.');
% quiver3(xyz(l,1), xyz(l,2), xyz(l,3), d(1), d(2), d(3), 0, 'r', 'maxheadsize', 1);
% plot3(p(:,1), p(:,2), p(:,3), 'go--');


axis image vis3d
view(3)
camproj('persp');
camlight(30,30);

addpath('/stck/bandrieu/Bureau/CYPRES/Regression3D/GPeyre/toolbox_graph/');
e = compute_edges(tri)';
le = sqrt(sum( (xyz(e(:,1),:) - xyz(e(:,2),:)).^2, 2));
figure;
hist(le)

mean(le)
min(le)
max(le)

