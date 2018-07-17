clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

cl = colorcet( 'I2', 'N', 2 );

uvxyz = [
    -1.0000000000000000       -1.0000000000000000       -1.0000000000000000        1.0000000000000000        1.0000000000000000      -0.11755705045849460       0.16180339887498874
    1.0000000000000000       -1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000       0.11755705045849460       0.16180339887498874
    ];


if 0
    
    n = 200;
    t = linspace( -1, 1, n )';
    
    
    figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
    hold on
    
    for isurf = 1:2
        b = readCoeffs2( sprintf('debugbsi_reg%d.bern', isurf) );
        
        s = eval_Bezier_patch_tpgrid(b, t, t);
        surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
            'edgecolor','none','facecolor',cl(isurf,:),'diffusestrength',0.7, ...
            'facelighting','gouraud','specularstrength',0.15);
    end
    
    axis image vis3d
    view(120,30)
    camproj('persp');
    camlight(30,30);
end

%%
uvbox = [
    -1.0000000000000000        1.0000000000000000       -1.0000000000000000        1.0000000000000000
    -1.0000000000000000        1.0000000000000000       -1.0000000000000000        1.0000000000000000
    ];

duv_ds = [
    0.50249378105604436       -4.5094136758017100E-015   5.0000000000000003E-002  -4.4867385072421875E-015
    0.50249378105604436        4.5094136758017100E-015   5.0000000000000003E-002   4.4867385072421875E-015
    ];



figure;
for isurf = 1:2
    subplot(1,2,isurf);
    hold on
    rectangle('position', [uvbox(isurf,[1,3]), uvbox(isurf,[2,4]) -  uvbox(isurf,[1,3])]);
    
    quiver(uvxyz(:,2*isurf-1), uvxyz(:,2*isurf), ...
        duv_ds(:,2*isurf-1), duv_ds(:,2*isurf), 'r');
    
    for j = 1:size(uvxyz,1)
        text(uvxyz(j,2*isurf-1), uvxyz(j,2*isurf), num2str(j), 'color', 'b');
    end
    
    axis image
end
