clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

cl = colorcet( 'I2', 'N', 2 );

uvxyz = [
    5.2506781322596369E-014  0.25833480542467829       -9.4987392425537109E-002  0.16535034245563257       -6.5981967386681408E-014   2.2945250472328091E-007  -2.4757973449140991E-014
   3.4433921973202482E-014  0.25833487510681152       -9.4987322743403949E-002  0.16535034245562349       -4.3270942394697298E-014   1.8566992919648613E-007  -1.5688204582205944E-014
    ];

for isurf = 1:2
        S(isurf).C = readCoeffs2( sprintf('debugssi_surf%d.cheb', isurf) );
        S(isurf).D = cheb_diff2( S(isurf).C );
        D2u = cheb_diff2( S(isurf).D(:,:,:,1) );
        D2v = cheb_diff2( S(isurf).D(:,:,:,2) );
        S(isurf).D2 = cat(4, D2u, D2v(:,:,:,2));
end


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
    0.0000000000000000        1.1920928955078125E-007  0.25833487510681152       0.25833499431610107     
  -9.4987392425537109E-002  -9.4987273216247559E-002  0.16535019874572754       0.16535031795501709
    ];

% duv_ds = [
%     -17.110353517320299       -1.1880586312689108E-010  -17.110353517216026       -1.8424195438433798E-011
%     -15.527523788149484       -2.3174108536945138E-012  -15.527523788144286        2.0978743372413412E-012
% %     0.50249378105604436       -4.5094136758017100E-015   5.0000000000000003E-002  -4.4867385072421875E-015
% %     0.50249378105604436        4.5094136758017100E-015   5.0000000000000003E-002   4.4867385072421875E-015
%     ];

duv_ds = zeros(size(uvxyz,1),4);
for i = 1:size(uvxyz,1)
    i
    for isurf = 1:2
        isurf
        [ uvxyz(i,2*isurf-1:2*isurf) - uvbox(isurf,[1,3]), 
        uvbox(isurf,[2,4]) - uvxyz(i,2*isurf-1:2*isurf) ]
    end
    [ ~, uv_s ] = diffgeom_intersection_curve(S, reshape(uvxyz(i,1:4), 2, 2));
    duv_ds(i,:) = reshape(uv_s(:,1,:), 1, 4);
end



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
