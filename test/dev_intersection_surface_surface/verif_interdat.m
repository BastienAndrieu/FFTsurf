clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

fid = fopen( 'interdat_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );

for ic = 1:nc
    endpoints = str2num( fgetl( fid ) );
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


figure;
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:) );
end

for ic = 1:nc
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), 'y.-')
end

axis image vis3d
view(3)
camlight(30,30);
camproj('persp');