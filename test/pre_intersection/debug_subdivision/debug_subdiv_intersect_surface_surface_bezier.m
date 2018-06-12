clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

cl = colorcet( 'I2', 'N', 10 );
n = 200;
t = linspace(0,1,n)';


for isurf = 2
    tree = importdata(sprintf('surf_%d_tree.dat',isurf));
    umin = min(tree(:,1));
    vmin = min(tree(:,3));
    umax = max(tree(:,2));
    vmax = max(tree(:,4));
    
    root.c = readCoeffs2(sprintf('surf_%d_root.cheb',isurf));
    parent.b = readCoeffs2(sprintf('surf_%d_parent.bern',isurf));
    parent.c = chgvar2( root.c, [umin, vmin ; umax, vmax] );
    
    for ichild = 1:4
        child(ichild).b = readCoeffs2(sprintf('surf_%d_child_%d.bern',isurf,ichild));
        child(ichild).c = chgvar2( root.c, [tree(ichild,[1,3]) ; tree(ichild,[2,4])] );
    end
    
    
    
    figure;
    hold on
    
    for j = 1:size(tree,1)
        rectangle( 'position', ...
            [ tree(j,[1,3]), tree(j,[2,4]) - tree(j,[1,3]) ] );
        
        %         xm = 0.5*sum(tree(j,1:2));
        %         ym = 0.5*sum(tree(j,3:4));
        %         text( xm, ym, num2str(j), 'fontsize', 6, 'color', 'r' );
    end
%     rp = rectangle( 'position', [umin, vmin, umax-umin, vmax-vmin] );
%     set( rp, 'edgecolor', 'r', 'linewidth', 2, 'linestyle', ':' );
    
    axis( repmat([-1,1],1,2) );
    daspect( [1,1,1] );
    
    
    
    figure;
    hold on
    
    [sroot, lu, lv] = surf_chebyshev( root.c, 1, 200, 10, 0, 0 );
    set( sroot, 'facecolor', 'none' );
    set( [lu;lv], 'color', 0.5*[1,1,1] );
    
    %     [sparent,~,~,~,lc] = surf_chebyshev( parent.c, 1, 200, 0, 0, 1 );
    %     set( sparent, 'facecolor', 'none' );
    s = eval_Bezier_patch_tpgrid( parent.b, t, t );
    lc = contourSurf(s(:,:,1),s(:,:,2),s(:,:,3));
    set( lc, 'linewidth', 1.5, 'linestyle', '--' );
    
    for ichild = 1:4
        s = eval_Bezier_patch_tpgrid( child(ichild).b, t, t );
        surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
            'edgecolor','none','facecolor',cl((isurf-1)*2 + ichild, :),'diffusestrength',0.7, ...
            'facelighting','gouraud','specularstrength',0.15);
        %         schild = surf_chebyshev( child(ichild).c, 1 );
        %         set( schild, 'facecolor', cl((isurf-1)*2 + ichild, :) );
    end
    
    axis image vis3d
    view(3)
    camlight(30,30);
    camproj('persp');
    
    
    
    
    figure;
    hold on
    
    [sroot, lu, lv] = surf_chebyshev( root.c, 1, 200, 10, 0, 0 );
    set( sroot, 'facecolor', 'none' );
    set( [lu;lv], 'color', 0.5*[1,1,1] );
    
    [sparent,~,~,~,lc] = surf_chebyshev( parent.c, 1, 200, 0, 0, 1 );
    set( sparent, 'facecolor', 'none' );
    set( lc, 'linewidth', 1.5, 'linestyle', '--' );
    
    for ichild = 1:4
        schild = surf_chebyshev( child(ichild).c, 1 );
        set( schild, 'facecolor', cl((isurf-1)*2 + ichild, :) );
    end
    
    axis image vis3d
    view(3)
    camlight(30,30);
    camproj('persp');
    
    
end