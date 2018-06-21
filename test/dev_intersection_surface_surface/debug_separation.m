clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');

vec = importdata('vec.dat')';

R = 0;
for isurf = 1:2
    pn(isurf).b = readCoeffs2(sprintf('bpn%d.bern',isurf));
    pn(isurf).sep = importdata(sprintf('sep%d.dat',isurf));
    R = max( R, max(max( sum(pn(isurf).b.^2,3 ) ) ) );
    
    pn(isurf).s = pn(isurf).b ./ repmat( sqrt(sum( pn(isurf).b.^2, 3 )),[1,1,3] );
end
% R = sqrt(R);
R = 1.1;

cl = ['r','b'];

figure;
hold on

for isurf = 1:2
    if length( pn(isurf).b(:,:,1) ) > 3
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

plot_separating_plane( vec );

axis( repmat( R*[-1,1], 1, 2 ) );
daspect([1,1,1]);
axis vis3d
view(3)
camlight(30,30);

%%
figure;
hold on

plot( 0, 0, 'kx' );

for isurf = 1:2
    if ~isempty( pn(isurf).sep )
        plot( pn(isurf).sep(:,2),  pn(isurf).sep(:,3), '.', 'color', cl(isurf) );
    end
end

axis image


