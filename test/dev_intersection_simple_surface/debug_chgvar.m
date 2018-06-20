clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

rootc = readCoeffs2( 'root_s.cheb' );
uvbox = [-1.0000000000000000      -0.50000000000000000        0.0000000000000000       0.50000000000000000];

regionc = chgvar2( rootc, reshape(uvbox,2,2) );

c = readCoeffs2( 'region_s.cheb' );
b = readCoeffs2( 'region_s.bern' );

figure;
hold on


surf_chebyshev( c, 1 );

surf( ifcht2d( regionc(:,:,1) ), ifcht2d( regionc(:,:,2) ), ifcht2d( regionc(:,:,3) ), ...
    'facecolor', 'none', 'edgecolor', 'b' );

surf( b(:,:,1), b(:,:,2), b(:,:,3), 'facecolor', 'none', 'edgecolor', 'r' );

axis image vis3d
% axis( xm([1,1,2,2,3,3]) + 10*r*repmat([-1,1],1,3) );
daspect( [1,1,1] );
view(3)
camproj('persp');
camlight(30,30);

