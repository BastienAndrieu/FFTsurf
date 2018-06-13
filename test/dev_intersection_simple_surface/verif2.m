clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

for isurf = 1:2
    surfa(isurf).c = readCoeffs2( sprintf('c_%d.cheb',isurf) );
end

uvxyz = importdata('uv_xyz.dat');

cls = colorcet( 'I2', 'N', 2 );

figure;
hold on

for isurf = 1:2
    si = surf_chebyshev( surfa(isurf).c, 1, 300 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
end

if ~isempty(uvxyz)
    plot3( uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), 'r*', 'markersize', 4 );
end

axis image vis3d
view(3)
% view(-130,28)
camproj('persp');
camlight(30,30);




% 
% b1 = chebyshev2bezier_2( surfa(2).c );
% b2 = readCoeffs2('root_s_bezier.bern');
% 
% b1 - b2

return

%%
n = 500;
t = linspace( 0, 1, n )';
u = linspace( -1, 1, n )';

c = surfa(1).c;
b = readCoeffs2('root_s_bezier.bern');

sc = ICT2tpgrid( c, u, u );
sb = eval_Bezier_patch_tpgrid( b, t, t );


e = sc - sb;

figure;
surf( log10( sqrt(sum(e.^2,3)) ), ...
    'facecolor', 'interp', 'edgecolor', 'none' );
view(2);
axis image
colorbar