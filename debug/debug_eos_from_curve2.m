clc; clear; close all
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');



v = linspace(0,1,10);

for i = 1:2
    erl{i} = importdata(sprintf('eos_from_curve2_erl_%d.dat', i));
end

xyz = importdata('eos_from_curve2_xyz.dat');
dxyz_dw = importdata('eos_from_curve2_dxyz_dw.dat');

center = importdata('eos_from_polyline_center.dat');
localframe = importdata('eos_from_polyline_localframe.dat');
axes = reshape(localframe(:,1:6), [size(localframe,1),3,2]);
angle = localframe(:,7);

s = readCoeffs2('eos_from_polyline_cgl_tpgrid.dat');


figure;
hold on

plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'k.-');
quiver3(xyz(:,1), xyz(:,2), xyz(:,3), ...
    dxyz_dw(:,1), dxyz_dw(:,2), dxyz_dw(:,3), 'r');

for i = 1:2
    plot3(erl{i}(:,1), erl{i}(:,2), erl{i}(:,3), 'o-');
end


for i = 1%:2
    quiver3(center(:,1), center(:,2), center(:,3), ...
        axes(:,1,i), axes(:,2,i), axes(:,3,i), 0);
end

for j = 1:numel(angle);
    a = angle(j)*v;
    plot3(...
        center(j,1) + cos(a)*axes(j,1,1) + sin(a)*axes(j,1,2), ...
        center(j,2) + cos(a)*axes(j,2,1) + sin(a)*axes(j,2,2), ...
        center(j,3) + cos(a)*axes(j,3,1) + sin(a)*axes(j,3,2), ...
        'b-');
end

surf(s(:,:,1), s(:,:,2), s(:,:,3), ...
    'facecolor', 'none', 'facelighting', 'gouraud', 'specularstrength', 0);

axis image vis3d
view(130,15);
camproj('persp');
camlight(30,30);


% w = importdata('eos_from_curve2_w.dat');
% figure;
% plot(w,'o-');
