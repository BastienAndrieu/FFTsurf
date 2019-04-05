clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

c = readCoeffs2('c.cheb');

uv = importdata('geodesic_uv.dat');

xyz = ICT2unstr(c, uv);


%%
niso = 9;
lwiso = 0.1;

figure('units','normalized', 'position', [.1 .1 .8 .8]);

subplot(1,2,1)
hold on

[S, lu, lv] = surf_chebyshev(c, 1, 200, 9);
set([lu;lv], 'linewidth', lwiso);

plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'k.-', 'linewidth', 1, 'markersize', 5)

axis image vis3d
view(23,16)
camlight(30,30);
camproj('persp');


subplot(1,2,2);
hold on

tiso = linspace(-1,1,niso);
for i = 1:niso
    line(tiso(i)*[1,1], [-1,1], 'color', 'r', 'linewidth', lwiso);
    line([-1,1], tiso(i)*[1,1], 'color', 'b', 'linewidth', lwiso);
end

plot(uv(:,1), uv(:,2), 'k.-', 'linewidth', 1, 'markersize', 5)

axis(repmat([-1,1], 1, 2));
daspect([1,1,1]);
