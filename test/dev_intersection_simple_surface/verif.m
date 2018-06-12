clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');


bc = readCoeffs1('root_c_bezier.bern');
bs = readCoeffs2('root_s_bezier.bern');

tuvxyz = importdata('tuv_xyz.dat');
treec = importdata('tree_c.dat');
trees = importdata('tree_s.dat');


n = 100;
t = linspace( 0, 1, n )';



figure;
hold on

s = eval_Bezier_patch_tpgrid( bs, t, t );
surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
    'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
    'facelighting','gouraud','specularstrength',0.15);
% surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), ...
%     'facecolor', 'none', 'edgecolor', 'r', 'linestyle', '--' );


c = evalBernstein( bc, t );
plot3( c(:,1), c(:,2), c(:,3), 'k-' );
% plot3( bc(:,1), bc(:,2), bc(:,3), 'b.--' );

if ~isempty(tuvxyz)
    plot3( tuvxyz(:,4), tuvxyz(:,5), tuvxyz(:,6), 'r*', 'markersize', 4 );
end

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);


%%
figure;

subplot(1,2,1);
hold on
for i = 1:size(treec,1)
    plot( treec(i,:), [0,0], 'k+-' );
end
if ~isempty(tuvxyz)
    plot( tuvxyz(:,1), 0*tuvxyz(:,1), 'r*', 'markersize', 4 );
end
axis([-1,1,-0.2,0.2]);
daspect([1,1,1]);

subplot(1,2,2);
hold on
for i = 1:size(trees,1)
    rectangle( 'position', [trees(i,1), trees(i,3), ...
        trees(i,2) - trees(i,1), trees(i,4) - trees(i,3) ] );
end
if ~isempty(tuvxyz)
    plot( tuvxyz(:,2), tuvxyz(:,3), 'r*', 'markersize', 4 );
end
axis image