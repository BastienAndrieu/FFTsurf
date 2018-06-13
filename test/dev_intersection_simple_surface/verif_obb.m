clc; clear; close all

cl = ['b','r'];

bc = readCoeffs1('region_c_bezier.bern');
bs = readCoeffs2('region_s_bezier.bern');



n = 300;
t = linspace( 0, 1, n )';

figure;
hold on

s = eval_Bezier_patch_tpgrid( bs, t, t );
surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
    'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
    'facelighting','gouraud','specularstrength',0.15);

c = evalBernstein( bc, t );
plot3( c(:,1), c(:,2), c(:,3), 'k-' );

for i = 1:2
    
    if i == 1
        box = importdata('xyzbox_c.dat');
    else
        box = importdata('xyzbox_s.dat');
    end
    ctr = box(1,:)';
    rng = box(2,:)';
    axe = box(3:5,:)';
    
    h = plotOBB( ctr, rng, axe );
    set( h, 'facecolor', 'none', 'edgecolor', cl(i) );
end

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);