clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');


instant = 10;

tri = importdata('../connect_01.dat');
xkp1 = importdata(sprintf('../pos_%3.3d.dat',instant+1));
xk = importdata(sprintf('../pos_%3.3d.dat',instant));
xkm1 = importdata(sprintf('../pos_%3.3d.dat',instant-1));

acc = xkp1 - 2*xk + xkm1;
mag_acc = sqrt(sum(acc.^2,2));

figure('units','normalized', 'position', [.1 .1 .8 .8]);
hold on

trisurf(tri, xk(:,1), xk(:,2), xk(:,3), mag_acc, 'facecolor', 'interp');

axis image vis3d
view(3)
% view(240,-50)
camproj('persp');





%%
xkp1 = importdata(sprintf('../pos_%3.3d.dat',instant+2));
xk = importdata(sprintf('../pos_%3.3d.dat',instant));
xkm1 = importdata(sprintf('../pos_%3.3d.dat',instant-2));

acc = xkp1 - 2*xk + xkm1;
mag_acc = 0.25 * sqrt(sum(acc.^2,2));

figure('units','normalized', 'position', [.1 .1 .8 .8]);
hold on

trisurf(tri, xk(:,1), xk(:,2), xk(:,3), mag_acc, 'facecolor', 'interp');

axis image vis3d
view(3)
% view(240,-50)
camproj('persp');


















