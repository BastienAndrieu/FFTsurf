clc; clear; close all


addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');


cl = colorcet( 'I2', 'N', 2 );
n = 200;
t = linspace(0,1,n)';



figure;
hold on


for isurf = 1:2
    b = readCoeffs2(sprintf('bs_debugobb_%d.bern',isurf));
    surf( b(:,:,1), b(:,:,2), b(:,:,3), ...
        'edgecolor',0.5*cl(isurf,:),'facecolor','none','linewidth',0.2);
    
    s = eval_Bezier_patch_tpgrid( b, t, t );
    surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
        'edgecolor','none','facecolor',cl(isurf,:),'diffusestrength',0.7, ...
        'facelighting','gouraud','specularstrength',0.15);
    
    
    if 1
        box = importdata(sprintf('xyzbox_debugobb_%d.dat',isurf));
        ctr = box(1,:)';
        rng = box(2,:)';
        axe = box(3:5,:)';
    else
        box = OBB_Bernstein2( b );
        ctr = box.center;
        rng = box.extents;
        axe = box.axes;
    end
    
    h = plotOBB( ctr, rng, axe );
    set( h, 'facecolor', 'none', 'edgecolor', cl(isurf,:) );
    
    
end


axis image vis3d
view(3)
camlight(30,30);
camproj('persp');



%%

clear box
for isurf = 1:2
    obb = importdata(sprintf('xyzbox_debugobb_%d.dat',isurf));
    box(isurf).center = obb(1,:)';
    box(isurf).extents = obb(2,:)';
    box(isurf).axes = obb(3:5,:)';
end

% box(2).center = box(2).center + 2*(2*rand(3,1)-1);



intersect_OBBs( box(1), box(2) )



figure;
hold on

for isurf = 1:2
    h = plotOBB( box(isurf).center, box(isurf).extents, box(isurf).axes );
    set( h, 'facecolor', cl(isurf,:), 'edgecolor', 'k' );
end

axis image vis3d
view(3)
camlight(30,30);
camproj('persp');
