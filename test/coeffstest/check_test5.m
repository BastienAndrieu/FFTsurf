clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
format long

% verification de la coincidences des coins des deux surfaces du test 5
% (cone/cylindre)

u = [-1.0, -1.0;
    1.0, -1.0;
    -1.0, 1.0;
    1.0, 1.0];

for isurf = 1:2
    S(isurf).c = readCoeffs2( sprintf( 'C%d_test05.txt', isurf ) );
%     S(isurf).c = readCoeffs2( ...
%         sprintf('/d/bandrieu/GitHub/FFTsurf/test/dev_intersection_surface_surface/c_%d.cheb',isurf) );
    S(isurf).b = chebyshev2bezier_2( S(isurf).c );
    S(isurf).x = ICT2unstr( S(isurf).c, u );
    [m,n,~] = size( S(isurf).b );
    for k = 1:2
        for j = 1:2
            S(isurf).y((k-1)*2+j,:) = reshape( ...
                S(isurf).b( (j-1)*(m-1) + 1, (k-1)*(n-1) + 1, 1:3 ), 1, 3 );
        end
    end
    
%     S(isurf).x - S(isurf).y
end

for j = 1:4
    for i = 1:4
        fprintf('i = %d, j = %d, dc = %e\t db = %e\n', ...
            i, j, ...
            norm( S(1).x(i,:) - S(2).x(j,:) ), ...
            norm( S(1).y(i,:) - S(2).y(j,:) ) );
    end
end



cl = colorcet( 'I2', 'N', 2 );
cls = CC( cl, 0.0, 0.8, 2.0 );
clb = CC( cl, 0.0, 1.5, 0.5 );

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on
for isurf = 1:2
    si = surf_chebyshev( S(isurf).c, 1 );
    set( si, 'facecolor', cls(isurf,:), 'specularstrength', 0 );
    
%     surf( S(isurf).b(:,:,1), S(isurf).b(:,:,2), S(isurf).b(:,:,3), ...
%     'facecolor', 'none', 'edgecolor', clb(isurf,:), 'linestyle', '-' );


    quiver3( S(isurf).b(1,1,1), S(isurf).b(1,1,2), S(isurf).b(1,1,3), ...
        S(isurf).b(2,1,1) - S(isurf).b(1,1,1), ...
        S(isurf).b(2,1,2) - S(isurf).b(1,1,2), ...
        S(isurf).b(2,1,3) - S(isurf).b(1,1,3), 2, 'r' );
    quiver3( S(isurf).b(1,1,1), S(isurf).b(1,1,2), S(isurf).b(1,1,3), ...
        S(isurf).b(1,2,1) - S(isurf).b(1,1,1), ...
        S(isurf).b(1,2,2) - S(isurf).b(1,1,2), ...
        S(isurf).b(1,2,3) - S(isurf).b(1,1,3), 0.5, 'g' );
end

l = [1];
for isurf = 1:2
    plot3( S(isurf).x(l,1), S(isurf).x(l,2), S(isurf).x(l,3), '*', 'color', clb(isurf,:) );
end

sqrt( sum( ( S(1).x(l,:) - S(2).x(l,:) ).^2, 2 ) )


axis image vis3d
view(3)
camproj('persp');
camlight(30,30);