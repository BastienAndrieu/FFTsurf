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

plot3( points(:,1), points(:,2), points(:,3), 'k.', 'markersize', 8 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );



%%
return
uvxyz = [0.0000000000000000      -0.56311007890829679        6.7229545007077693E-014 -0.55511007890828157       0.28155503945414539       0.14142135623731170       0.14142135623730576 ];
uv = reshape( uvxyz(1:4), 2, 2 );
xyz = uvxyz(5:7);

plot3( xyz(1), xyz(2), xyz(3), 'ro', 'markersize', 8 );


for isurf = 1:2
    surfa(isurf).C = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    surfa(isurf).D = cheb_diff2( surfa(isurf).C );
    surfa(isurf).D2u = cheb_diff2( surfa(isurf).D(:,:,:,1) );
    surfa(isurf).D2v = cheb_diff2( surfa(isurf).D(:,:,:,2) );
    
    %     xyzi = ICT2unstr( surfa(isurf).C, uv(:,isurf)' );
    %     plot3( xyzi(:,1), xyzi(:,2), xyzi(:,3), 'k*' );
end

[dxyz_ds, duv_ds] = geomdiff_intersection_curve( surfa, uv )


quiver3( xyz(1), xyz(2), xyz(3), dxyz_ds(1), dxyz_ds(2), dxyz_ds(3), 'k' );


uvbox = [...
    -1.0000000000000000        0.0000000000000000       -1.0000000000000000      -0.50400000000009904
    0.0000000000000000        1.0000000000000000       -1.0000000000000000      -0.50000000000009837];

figure;
for isurf = 1:2
    subplot(1,2,isurf);
    hold on
    rectangle( 'position', ...
        [ uvbox(isurf,[1,3]), ...
        uvbox(isurf,[2,4]) - uvbox(isurf,[1,3]) ] );
    
    quiver( uv(1,isurf), uv(2,isurf), duv_ds(1,isurf), duv_ds(2,isurf), 0.1 );
    axis( repmat( 1.01*[-1,1], 1, 3 ) );
    daspect( [1,1,1] );
end
