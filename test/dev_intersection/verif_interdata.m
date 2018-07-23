clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

nsurf = 2;
PLOT_TREE = 1;

cl = colorcet( 'I2', 'N', nsurf );
cl = CC( cl, 0.0, 0.8, 1.5 );
cli = 'k';%'y';

%% read interdata

points = importdata( 'interdataglobal_points.dat' );

fid = fopen( 'interdataglobal_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );

for ic = 1:nc
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
end
fclose( fid );

%% Plot

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for isurf = 1:nsurf
    c = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    
    if PLOT_TREE
        tree = importdata(sprintf('treessi_%d.dat',isurf));
        for i = 1:size(tree,1)
            m = round( max(tree(i,[2,4]) - tree(i,[1,3])) * 50 );
            uvi = [ linspace(tree(i,1),tree(i,2),m)' , tree(i,3)*ones(m,1) ;
                tree(i,2)*ones(m,1), linspace(tree(i,3),tree(i,4),m)' ;
                linspace(tree(i,2),tree(i,1),m)' , tree(i,4)*ones(m,1) ;
                tree(i,1)*ones(m,1), linspace(tree(i,4),tree(i,3),m)' ];
            
            l = ICT2unstr( c, uvi );
            
            plot3( l(:,1), l(:,2), l(:,3), '-', ...
                'color', 0.4*cl(isurf,:), 'linewidth', 0.1 );
        end
    end
end

for ic = 1:nc
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
        '-', 'color', cli, 'linewidth', 1.0, 'markersize', 5 )
end

plot3( points(:,1), points(:,2), points(:,3), '.', 'color', cli, 'markersize', 10 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );



%%

ncs = 0;
for ic = 1:nc
    if curves(ic).nsplit > 0
        ncs = ncs + curves(ic).nsplit - 1;
    end
end

clrc = hsv(ncs);
clrc = clrc(randperm(ncs),:);

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

plot3( points(:,1), points(:,2), points(:,3), '.', 'color', cli, 'markersize', 10 );

js = 0;
for ic = 1:nc
    if curves(ic).nsplit > 0
        for is = 1:curves(ic).nsplit - 1
            js = js + 1;
            l = curves(ic).isplit(is,2):curves(ic).isplit(is+1,2);
            plot3( curves(ic).xyz(l,1), curves(ic).xyz(l,2), curves(ic).xyz(l,3), ...
                '-', 'color', clrc(js,:), 'markersize', 5 )
        end
    else
        plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
            '-', 'color', cli, 'markersize', 5 )
    end
end

axis image vis3d
view(3)
camproj('persp');

%%
figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for ic = 1:nc
    plot3( curves(ic).xyz(:,1), curves(ic).xyz(:,2), curves(ic).xyz(:,3), ...
        '.-', 'color', cli, 'markersize', 5 )
end

axis image vis3d
view(3)
camproj('persp');





