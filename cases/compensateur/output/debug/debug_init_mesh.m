clc; clear; close all

tri = importdata('../connect_01.dat');
xyz = importdata('../pos_000.dat');

% tri = importdata('../../init/mesh/trib.dat');
% xyz = importdata('../../init/mesh/xyzb.dat');

fprintf('%d verts, %d triangles\n',size(xyz,1), size(tri,1));

figure;
hold on
trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'y', 'specularstrength', 0);
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

return

x = [
    5.00E-03
    0.00031618
    0.00028972
    2.8386e-04
    2.00E-04
    1.00E-04
    ];
y = [
    2826
5082
5570
5681
8160
16294
];

p = polyfit(1./x, y, 2);
q = p;
q(3) = q(3) - 23238;

r = roots(q);
r(r < 0) = [];
1./r
