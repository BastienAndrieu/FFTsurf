clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');


for isurf = 1:2
    S(isurf).x  = readCoeffs2( sprintf( 'region_%d_x.bern', isurf ) );
    S(isurf).pn = readCoeffs2( sprintf( 'region_%d_pn.bern', isurf ) );
    
    S(isurf).sep = reshape( S(isurf).pn, [], 3 );
end

cl = linspecerp(2);
n = 200;
t = linspace(0,1,n)';


%% Surfaces
figure;
hold on

for isurf = 1:2
    s = eval_Bezier_patch_tpgrid( S(isurf).x, t, t );
    surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
        'edgecolor', 'none', 'facecolor',cl(isurf, :), 'diffusestrength',0.7, ...
        'facelighting', 'gouraud', 'specularstrength', 0.15);
    
    surf( S(isurf).x(:,:,1), S(isurf).x(:,:,2), S(isurf).x(:,:,3), ...
        'facecolor', 'none', 'edgecolor', 0.5*cl(isurf,:) )
end

axis image vis3d
view(3)
camlight(30,30);
camproj('persp');



%% Gauss maps
figure;
hold on

for isurf = 1:2
    if length(S(isurf).pn) == 3
        quiver3( 0, 0, 0, S(isurf).pn(1), S(isurf).pn(2), S(isurf).pn(3), 0, 'color', 0.5*cl(isurf,:) );
    else
        s = eval_Bezier_patch_tpgrid( S(isurf).pn, t, t );
        surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
            'edgecolor', 'none', 'facecolor',cl(isurf, :), 'diffusestrength',0.7, ...
            'facelighting', 'gouraud', 'specularstrength', 0.15);
        
        surf( S(isurf).pn(:,:,1), S(isurf).pn(:,:,2), S(isurf).pn(:,:,3), ...
            'facecolor', 'none', 'edgecolor', 0.5*cl(isurf,:) )
    end
end

axis image vis3d
view(3)
camlight(30,30);
camproj('persp');


%% Gauss map separation
[ p1, stat1 ] = separating_plane( S(1).sep, S(2).sep )
[ p2, stat2 ] = separating_plane( S(1).sep, -S(2).sep )

figure;
hold on

for isurf = 1:2
%     bpn = [ [0,0,0]; S(isurf).sep ];
%     [ ep, fp ] = plot_spherical_pyramid( ...
%     bpn, convhulln(bpn), 1, 1.0 );
%     set( fp, 'facecolor', cls(isurf,:) );
%     set( ep, 'color', 'k', 'linewidth', 1 );
%     
%     sep = importdata( sprintf('sep%d.dat',isurf) );
%     size(sep)
%     
%     sep = sep ./ repmat( sqrt( sum( sep.^2, 2 ) ), 1, 3 );
%     %     plot3( sep(:,1), sep(:,2), sep(:,3), '.', 'color', 0.7*cls(isurf,:) );
%     sep = reshape( sep, [size(S(isurf).pn,1), size(S(isurf).pn,2), 3] );
%     surf( sep(:,:,1), sep(:,:,2), sep(:,:,3), ...
%         'facecolor', 'none', 'edgecolor', 0.5*cls(isurf,:) );
    sep = S(isurf).sep ./ repmat( sqrt( sum( S(isurf).sep.^2, 2 ) ), 1, 3 );
    plot3( sep(:,1), sep(:,2), sep(:,3), '.', 'color', 0.7*cl(isurf,:) );
end

if stat1 == 0
[ q, p ] = plot_separating_plane( p1 );
end

if stat2 == 0
[ q, p ] = plot_separating_plane( p2 );
set( p, 'facecolor', 'g' );
end

axis image vis3d
view(3)
camlight(30,30);