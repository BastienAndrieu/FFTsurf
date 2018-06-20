clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

format long

cl = ['b','r'];


ccroot = readCoeffs1( 'root_c_x.cheb' );
% csroot = readCoeffs2( 'root_s_x.cheb' );
csroot = readCoeffs2( 'c_2.cheb' );

% tbox = [0.92397689819335938       0.92397880554199219];
% uvbox = [-0.51086044311523438      -0.51085853576660156       0.13840103149414062       0.13840293884277344];
tbox = [0.92382812500000000       0.92578125000000000];
uvbox = [-0.51171875000000000      -0.50976562500000000       0.13671875000000000       0.13867187500000000];

tuv = [0.92480468750000000      -0.51074218750000000       0.13769531250000000];

% cc = chgvar1( ccroot, tbox );
cs = chgvar2( csroot, reshape(uvbox,2,2) );
cc = ccroot;
% cs = csroot;

xc = chebval1( ccroot, tuv(1) );
xs = ICT2unstr( csroot, tuv(2:3) );

bc = readCoeffs1( 'region_c_bezier.bern' );
bs = readCoeffs2( 'region_s_bezier.bern' );
bsroot = readCoeffs2( 'root_s_bezier.bern' );
bcroot = readCoeffs1( 'root_c_bezier.bern' );

x = [ bc; reshape(bs,[],3) ];
nx = size(x,1);
xm = sum( x, 1 ) / nx;
r = sqrt( max( sum( ( x - repmat(xm,nx,1) ).^2, 2 ) ) );


n = 300;
t = linspace( -1, 1, n )';

cls = colorcet( 'I2', 'N', 2 );

figure;
hold on

% for isurf = 1:2
%     c = readCoeffs2( sprintf('c_%d.cheb',isurf) );
% 
%     si = surf_chebyshev( c, 1, n );
%     set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
% end
surf_chebyshev( cs, 1, n );
surf( bsroot(:,:,1), bsroot(:,:,2), bsroot(:,:,3), 'facecolor', 'none', 'edgecolor', 'r' );

g = chebval1( cc, t );
plot3( g(:,1), g(:,2), g(:,3), 'k' );

% plot3( ifcht2d( cs(:,:,1) ), ifcht2d( cs(:,:,2) ), ifcht2d( cs(:,:,3) ), 'r.' );
% plot3( ifcht( cc(:,1) ), ifcht( cc(:,2) ), ifcht( cc(:,3) ), 'b.' );


% surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), 'facecolor', 'none', 'edgecolor', 'r' );
% plot3( bs(:,:,1), bs(:,:,2), bs(:,:,3), 'r*' );
% plot3( bc(:,1), bc(:,2), bc(:,3), 'b*' );
plot3( bcroot(:,1), bcroot(:,2), bcroot(:,3), 'b-' );


axis image vis3d
% axis( xm([1,1,2,2,3,3]) + 10*r*repmat([-1,1],1,3) );
daspect( [1,1,1] );
view(3)
camproj('persp');
camlight(30,30);

