clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

% pair = [7,37];
uv = [
    -0.11947839735255582       0.99999999999965250      
    -0.99999999999969180      -0.99999999999889233];

cl = colorcet( 'I2', 'N', 2 );
cl = CC( cl, 0.0, 0.8, 1.5 );


figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

x = zeros(2,3);
for isurf = 1:2
%     c = readCoeffs2( sprintf('../Jouke/propergol/C_%3.3d.cheb', pair(isurf)) );
    c = readCoeffs2( sprintf('surfroot%2.2d_x.cheb', isurf) );
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );    
 
    d = cheb_diff2(c);
    
    x(isurf,:) = ICT2unstr(c, uv(isurf,:));
    xu = ICT2unstr(d(:,:,:,1), uv(isurf,:));
    xv = ICT2unstr(d(:,:,:,2), uv(isurf,:));
    n = cross(xu,xv);
    n = n / norm(n);
    
    quiver3(x(isurf,1), x(isurf,2), x(isurf,3), n(1), n(2), n(3), 'color', 0.5*cl(isurf,:));
    
end

plot3(x(:,1), x(:,2), x(:,3), 'k.--');

axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );