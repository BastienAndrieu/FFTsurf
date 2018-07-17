clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

cl = colorcet( 'I2', 'N', 2 );

bc = readCoeffs1('debugbsi_c.bern');
bs = readCoeffs2('debugbsi_s.bern');

cc = readCoeffs1('debugbsi_c.cheb');
cs = readCoeffs2('debugbsi_s.cheb');

tuv = [
    -0.83333333333333326      -0.93333333333333335      -0.91666666666666663
    ];

n = 200;
t = linspace( -1, 1, n )';

figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

s = surf_chebyshev(cs, 1, n);
set(s, 'specularstrength', 0);

c = chebval1(cc, t);
plot3( c(:,1), c(:,2), c(:,3), 'k-' );

% s = eval_Bezier_patch_tpgrid( bs, t, t );
% surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
%     'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
%     'facelighting','gouraud','specularstrength',0.15);
%
% c = evalBernstein( bc, t );
% plot3( c(:,1), c(:,2), c(:,3), 'k-' );

s = icht2_pc(cs, tuv(:,2:3));
plot3( s(:,1), s(:,2), s(:,3), 'ro' );

c = chebval1(cc, tuv(:,1));
plot3( c(:,1), c(:,2), c(:,3), 'b.' );

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);