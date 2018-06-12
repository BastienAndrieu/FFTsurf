clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

cl = colorcet( 'I2', 'N', 10 );


for isurf = 2
    root.c = readCoeffs2(sprintf('surf_%d_root.cheb',isurf));
    parent.c = readCoeffs2(sprintf('surf_%d_parent.cheb',isurf));
    for ichild = 1:4
        child(ichild).c = readCoeffs2(sprintf('surf_%d_child_%d.cheb',isurf,ichild));
    end
    
    tree = importdata(sprintf('surf_%d_tree.dat',isurf));
    
    
    figure;
    hold on
    
    for j = 1:size(tree,1)
        rectangle( 'position', ...
            [ tree(j,[1,3]), tree(j,[2,4]) - tree(j,[1,3]) ] );
        
        %         xm = 0.5*sum(tree(j,1:2));
        %         ym = 0.5*sum(tree(j,3:4));
        %         text( xm, ym, num2str(j), 'fontsize', 6, 'color', 'r' );
    end
    
    axis( repmat([-1,1],1,2) );
    daspect( [1,1,1] );
    
    
    
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