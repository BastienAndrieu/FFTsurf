clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

points = importdata( 'interdat_points.dat' );

fid = fopen( 'interdat_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );

for ic = 1:nc
    curves(ic).endpoints = str2num( fgetl( fid ) );
    curves(ic).uvbox = zeros( 2, 4 );
    for isurf = 1:2
        curves(ic).uvbox(isurf,:) = str2num( fgetl( fid ) );
    end
    n = str2num( fgetl( fid ) );
    curves(ic).uv = zeros( n, 2, 2 );
    curves(ic).xyz = zeros( n, 3 );
    for ip = 1:n
        line = str2num( fgetl( fid ) );
        curves(ic).uv(ip,:,1) = line(1:2);
        curves(ic).uv(ip,:,2) = line(3:4);
        curves(ic).xyz(ip,:) = line(5:7);
    end
end
fclose( fid );



cl = colorcet( 'I2', 'N', 2 );
% cl = linspecerp(2);
cl = CC( cl, 0.0, 0.8, 1.5 );


figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
end

for ic = 1:nc
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), 'k-')
end

% plot3( points(:,1), points(:,2), points(:,3), 'k*' );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );





return


u = [0.0000000000000000      -0.50400000000009904;
    0.0000000000000000      -0.50000000000009837];
x = zeros(2,3);
for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    x(isurf,:) = ICT2unstr( c, u(isurf,:) );
    plot3( x(isurf,1), x(isurf,2), x(isurf,3), 'o', 'color', 0.15*cl(isurf,:) );
end

points(2,:)
x

