clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

fid = fopen( 'intersection_points.dat', 'r' );
np = str2num( fgetl(fid) );
xyz = zeros(np,3);
uv = zeros(np,4);
for i = 1:np
    line = str2num( fgetl(fid) );
    xyz(i,:) = line(1:3);
    uv(i,:) = line(4:7);
end
fclose(fid);


fid = fopen( 'intersection_curves.dat', 'r' );
nc = str2num( fgetl(fid) );
for i = 1:nc
    curv(i).endpoints = str2num( fgetl(fid) );
    curv(i).uvbox(1,:) = str2num( fgetl(fid) );
    curv(i).uvbox(2,:) = str2num( fgetl(fid) );
end
fclose(fid);



%% 3d
cls = colorcet( 'I2', 'N', 2 );
figure;
hold on


for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb',isurf) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:), 'facealpha', 0.2 );
end

plot3( xyz(:,1), xyz(:,2), xyz(:,3), 'k.', 'markersize', 8 )

for i = 1:nc
    plot3( xyz(curv(i).endpoints,1), ...
        xyz(curv(i).endpoints,2), ...
        xyz(curv(i).endpoints,3), ...
        'k-' );
end

axis image vis3d
daspect( [1,1,1] );
view(3)
camproj('persp');
camlight(30,30);




%% 2d
cl = linspecer(nc);

figure;

for isurf = 1:2
    subplot( 1, 2, isurf )
    hold on
    
    plot( uv(:,2*isurf-1), uv(:,2*isurf), 'k.', 'markersize', 5 )

    for i = 1:nc
        r = curv(i).uvbox(isurf,:);
        rectangle( 'position', [ r([1,3]), r([2,4]) - r([1,3]) ], ...
            'edgecolor', cl(i,:) );
        plot( uv(curv(i).endpoints,2*isurf-1), uv(curv(i).endpoints,2*isurf), ...
            '-', 'color', cl(i,:) );
    end
    axis image
end
