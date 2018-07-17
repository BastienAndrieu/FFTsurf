clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');


nsurf = 2;

cl = colorcet( 'I2', 'N', nsurf );
cl = CC( cl, 0.0, 0.8, 1.5 );
cli = 'k';%'y';%

PLOT_TREE = 1;





points = importdata( 'uvxyz.dat' );
xyz = points(:,5:7);

xi = zeros(size(points,1),3,2);



figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:nsurf
    c = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    
%     xi(:,:,isurf) = ICT2unstr(c, points(:,2*isurf-1:2*isurf));
    
    %     plot3( xi(:,1,isurf), xi(:,2,isurf), xi(:,3,isurf), 'o', 'color', 0.4*cl(isurf,:) );
    
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

plot3( xyz(:,1), xyz(:,2), xyz(:,3), '.', 'color', cli, 'markersize', 6 );


fid = fopen( 'interdata_curves.dat', 'r' );
nc = str2num( fgetl( fid ) );
curve = zeros(nc,2);
for ic = 1:nc
    endpoints = str2num( fgetl( fid ) );
    curve(ic,:) = endpoints;
    plot3( xyz(endpoints,1), xyz(endpoints,2), xyz(endpoints,3), ...
        '-', 'color', cli, 'linewidth', 1.25 );
    
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



% sqrt( sum( (xi(:,:,1) - xi(:,:,2)).^2, 2 ) )

%%
return


figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
for isurf = 1:2
    subplot(1,2,isurf)
    hold on
    
    uv = points(:,2*isurf-1:2*isurf);
    plot( uv(:,1), uv(:,2), 'k.', 'markersize', 6 );
    
    for ic = 1:nc
        plot( uv(curve(ic,:),1), uv(curve(ic,:),2), 'k-' );
    end
    
    axis( repmat([-1,1], 1, 2) );
    daspect([1,1,1]);
end



size(points,1) - nc
