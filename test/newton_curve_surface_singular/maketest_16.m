clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

EPS = 10.0 * eps('double');

format long

numtest = 16;
icurv = 1;
ivar = 1;
ival = 1;

S.c = readCoeffs2( sprintf( '../coeffstest/C%d_test%2.2d.txt', 1+mod(icurv,2), numtest ) );
S.c( abs(S.c) < EPS ) = 0;

S.d = cheb_diff2( S.c );
d2u = cheb_diff2( S.d(:,:,:,1) );
d2v = cheb_diff2( S.d(:,:,:,2) );
S.d2 = cat( 4, d2u, d2v(:,:,:,2) );

[Ms,Ns,~] = size( S.c );



c = readCoeffs2( sprintf( '../coeffstest/C%d_test%2.2d.txt', icurv, numtest ) );
G.c = border_parameterization( c, ivar, ival );
G.d = cheb_diff1( G.c );
G.d2 = cheb_diff1( G.d );


% G.c

% G.c(1,:) = G.c(1,:) + 0.1 * (2*rand(1,3) - 1);






%%
writeCoeffs2( S.c, 'surf.cheb' );
writeCoeffs1( G.c, 'curv.cheb' );


%% plots
n = 100;
t = linspace(-1,1,n)';


figure;
hold on

surf_chebyshev( S.c, 1 );

g = chebval1( G.c, t );

plot3( g(:,1), g(:,2), g(:,3), 'k-' );

axis image vis3d
view(3)
camproj( 'persp' );
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


%%
t = -8e-3;
uv = [1.0, 0.0];

g = chebval1( G.c, t );
s = ICT2unstr( S.c, uv );

norm( g - s )