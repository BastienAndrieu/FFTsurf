clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

format long

% bc = readCoeffs1('root_c_bezier.bern');
% bs = readCoeffs2('root_s_bezier.bern');

bc = readCoeffs1('region_c_bezier.bern');
bs = readCoeffs2('region_s_bezier.bern');

boxc = importdata('xyzbox_c.dat');
boxs = importdata('xyzbox_s.dat');


tuvxyz = importdata('tuv_xyz.dat');
% tuvxyz = [-0.90000000000000002        1.5005561735261401E-002  -1.0000000000000000        1.6672846372512762E-002
%  -0.90000000000000002      -0.78681730352931478       -1.0000000000000000      -0.87424144836590523];%
% tuvxyz = [1.66728459E-02   1.50055615E-02 -0.899999976       1.50055615E-02 -0.899999976       7.63869436E-13
%     0.857568622      0.771811724     -0.899999976      0.771811724     -0.899999976       1.56260005E-10
%     -0.874241471     -0.786817312     -0.899999976     -0.786817312     -0.899999976      -1.39424583E-11];


c = readCoeffs2('c_1.cheb');

treec = importdata('tree_c.dat');
trees = importdata('tree_s.dat');

m = 100;
n = 300;
t = linspace( 0, 1, n )';



figure;
hold on

s = eval_Bezier_patch_tpgrid( bs, t, t );
surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
    'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
    'facelighting','gouraud','specularstrength',0.15);
% surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), ...
%     'facecolor', 'none', 'edgecolor', 'r', 'linestyle', '--' );

if 1
    h = plotOBB( boxc(1,:)', boxc(2,:)', boxc(3:5,:)' );
    set( h, 'facecolor', 'none', 'edgecolor', 'b' );
    plot3( boxc(1,1), boxc(1,2), boxc(1,3), 'b*' );
    
    h = plotOBB( boxs(1,:)', boxs(2,:)', boxs(3:5,:)' );
    set( h, 'facecolor', 'none', 'edgecolor', 'r' );
end


if 0
    for i = 3%1:size(trees,1)
        0.5 * (trees(i,:) + 1)
        uvi = [ linspace(trees(i,1),trees(i,2),m)' , trees(i,3)*ones(m,1) ;
            trees(i,2)*ones(m,1), linspace(trees(i,3),trees(i,4),m)' ;
            linspace(trees(i,2),trees(i,1),m)' , trees(i,4)*ones(m,1) ;
            trees(i,1)*ones(m,1), linspace(trees(i,4),trees(i,3),m)' ];
        
        g = ICT2unstr( c, uvi );
        plot3( g(:,1), g(:,2), g(:,3), 'b-', 'linewidth', 0.2 );
    end
end


if 0
    for i = 1%:size(treec,1)
        tlim = 0.5 * (treec(i,:) + 1)
        c = evalBernstein( bc, linspace(tlim(1),tlim(2),m)' );
        plot3( c(:,1), c(:,2), c(:,3), 'k-' );
    end
else
    c = evalBernstein( bc, t );
    plot3( c(:,1), c(:,2), c(:,3), 'k-' );
    % plot3( bc(:,1), bc(:,2), bc(:,3), 'b.--' );
end

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




% return




%%
uv = [-0.90000000000000002      -0.78681730327869559       -1.0000000000000000      -0.87424144808743953
    -0.90000000000000002       0.77181174154253684       -1.0000000000000000       0.85756860171392990
    -0.90000000000000002        1.5005561737555132E-002  -1.0000000000000000        1.6672846375061168E-002];

test = 18;

cls = colorcet( 'I2', 'N', 2 );


xyz = zeros(size(uv,1),3,2);

figure;
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb',isurf) );
    xyz(:,:,isurf) = ICT2unstr( c, uv(:,2*isurf-[1,0]) );
    
    si = surf_chebyshev( c, 1 );
    set( si, 'specularstrength', 0, 'facecolor', cls(isurf,:) );
    
    %     plot3( xyz(:,1,isurf), xyz(:,2,isurf), xyz(:,3,isurf), 'k.', 'markersize', 6 );
end


c = evalBernstein( bc, t );
plot3( c(:,1), c(:,2), c(:,3), 'k-', 'linewidth', 1 );

axis image vis3d
%view(3)
view(-130,28)

camproj( 'persp' );
camlight(30,30);


% sqrt( sum( ( xyz(:,:,1) - xyz(:,:,2) ).^2, 2 ) )
