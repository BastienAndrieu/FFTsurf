clc; clear; close all

tri = importdata('../connect.dat');
xyz = importdata('../pos_000.dat');

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

