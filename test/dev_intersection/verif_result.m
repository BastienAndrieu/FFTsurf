clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

cl = colorcet( 'I2', 'N', 2 );
cl = CC( cl, 0.0, 0.8, 1.5 );

points = importdata( 'uvxyz.dat' );
xyz = points(:,5:7);


figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
end

plot3( xyz(:,1), xyz(:,2), xyz(:,3), 'k.', 'markersize', 6 );


fid = fopen( 'interdata_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );
for ic = 1:nc
    endpoints = str2num( fgetl( fid ) );
    plot3( xyz(endpoints,1), xyz(endpoints,2), xyz(endpoints,3), 'k-' );
    
    for isurf = 1:2
        str2num( fgetl( fid ) );
    end
    n = str2num( fgetl( fid ) );
    for ip = 1:n
        str2num( fgetl( fid ) );
    end
end
fclose( fid );

axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
