clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');

% vec = importdata('vec.dat')';
norcol = importdata('norcol.dat')';

R = 0;
for isurf = 1:2
    pn(isurf).b = readCoeffs2(sprintf('bpn%d.bern',isurf));
    
    if size(pn(isurf).b,1) == 1
        pn(isurf).b = repmat( pn(isurf).b, [2,1,1] );
    end
    if size(pn(isurf).b,2) == 1
        pn(isurf).b = repmat( pn(isurf).b, [1,2,1] );
    end
    
    pn(isurf).sep = importdata(sprintf('sep%d.dat',isurf));
    R = max( R, max(max( sum(pn(isurf).b.^2,3 ) ) ) );
    
    pn(isurf).s = pn(isurf).b ./ repmat( sqrt(sum( pn(isurf).b.^2, 3 )),[1,1,3] );
end
% R = sqrt(R);
R = 1.1;



 n = 200;
 t = linspace(0,1,n)';
 
 



cl = ['r','b'];
cls = 0.5 * ( ones(2,3) + [1,0,0;0,0,1] );

figure;
hold on

for isurf = 1:2
    if length( pn(isurf).b(:,:,1) ) > 3
%         s = eval_Bezier_patch_tpgrid( pn(isurf).b, t, t );
%         surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
%             'edgecolor','none','facecolor',cls(isurf,:),'diffusestrength',0.7, ...
%             'facelighting','gouraud','specularstrength',0.15);
        
        surf( pn(isurf).s(:,:,1), pn(isurf).s(:,:,2), pn(isurf).s(:,:,3), ...
            'facecolor', 'none', 'edgecolor', cl(isurf) );
%         surf( pn(isurf).b(:,:,1), pn(isurf).b(:,:,2), pn(isurf).b(:,:,3), ...
%             'facecolor', 'none', 'edgecolor', cl(isurf) );
%         plot3( pn(isurf).sep(:,1),  pn(isurf).sep(:,2),  pn(isurf).sep(:,3), ...
%             '.', 'color', cl(isurf) );
    else
        quiver3( 0, 0, 0, pn(isurf).b(1,1,1), pn(isurf).b(1,1,2), pn(isurf).b(1,1,3), ...
            0, 'color', cl(isurf) );
    end
end

% plot_separating_plane( vec );
quiver3( 0, 0, 0, norcol(1), norcol(1), norcol(3), 'k' );

axis( repmat( R*[-1,1], 1, 2 ) );
daspect([1,1,1]);
axis vis3d
view(3)
camlight(30,30);

% return

%%
[ vec, stat ] = separate_spherical_boxes( pn(1).sep, pn(2).sep, 1 )


figure;
hold on

plot( 0, 0, 'kx' );

for isurf = 1:2
    if ~isempty( pn(isurf).sep )
        plot( pn(isurf).sep(:,2),  pn(isurf).sep(:,3), '.', 'color', cl(isurf) );
    end
end


if stat == 0
    plot_separating_plane( vec([2,3,1]) );
end

axis image
view(2);

