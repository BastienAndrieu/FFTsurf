clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

nsurftotal = 8;
listsurf = 1:nsurftotal;
nsurf = length(listsurf);

cl = colorcet( 'I2', 'N', nsurftotal );
cl = cl(randperm(size(cl,1)),:);
cl = CC( cl, 0.0, 0.8, 1.5 );
cli = 'k';%'y';

%% read interdata

points = importdata( 'intersection_points.dat' );
np = size(points,1);

fid = fopen( 'intersection_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );


aabb = 1e6*repmat([1,-1],1,3);
for ic = 1:nc
    curves(ic).dummy = str2num( fgetl( fid ) );
    curves(ic).smooth = str2num( fgetl( fid ) );
    %curves(ic).endpoints = str2num( fgetl( fid ) );
    curves(ic).uvbox = zeros( 2, 4 );
    for isurf = 1:2
        curves(ic).uvbox(isurf,:) = str2num( fgetl( fid ) );
    end
    curves(ic).nsplit = str2num( fgetl( fid ) );
    curves(ic).isplit = zeros( curves(ic).nsplit, 2 );
    for ip = 1:curves(ic).nsplit
        curves(ic).isplit(ip,:) = str2num( fgetl( fid ) );
    end
    curves(ic).class = zeros( curves(ic).nsplit-1, 1 );
    for is = 1:curves(ic).nsplit-1
        curves(ic).class(is) = str2num( fgetl( fid ) );
    end
    n = str2num( fgetl( fid ) );
    curves(ic).np = n;
    curves(ic).uv = zeros( n, 2, 2 );
    curves(ic).xyz = zeros( n, 3 );
    for ip = 1:n
        line = str2num( fgetl( fid ) );
        curves(ic).uv(ip,:,1) = line(1:2);
        curves(ic).uv(ip,:,2) = line(3:4);
        curves(ic).xyz(ip,:) = line(5:7);
    end
    
    if n > 0%curves(ic).dummy
        aabb(1:2:5) = min(aabb(1:2:5), min(curves(ic).xyz,[],1));
        aabb(2:2:6) = max(aabb(2:2:6), max(curves(ic).xyz,[],1));
    end
end
fclose( fid );

scn_diameter = norm(aabb(2:2:6) - aabb(1:2:5));



%% Plot

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

% for isurf = listsurf
%     c = readCoeffs2( sprintf('../coef/c_%3.3d.cheb', isurf) );
%     
%     si = surf_chebyshev( c, 1, 100 );
%     set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
%     if isurf > 112
%         set( si, 'facecolor', 0.5*[1,1,1], 'facealpha', 0.4 );
%     end
% 
% 
% end

for ic = 1:nc
    %if curves(ic).dummy; continue; end
    if curves(ic).smooth; continue; end
    if isempty(curves(ic).xyz); continue; end
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
        '-', 'color', cli, 'linewidth', 1.0, 'markersize', 5 )
    e = curves(ic).isplit(:,1)';
    plot3( points(e,1), points(e,2), points(e,3), '.', 'color', cli, 'markersize', 10 );
end

% plot3( points(:,1), points(:,2), points(:,3), '.', 'color', cli, 'markersize', 10 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );