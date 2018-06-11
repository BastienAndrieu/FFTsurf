clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');


c = readCoeffs2('c.cheb');
b = readCoeffs2('b.bern');
child(1).b = readCoeffs2('bsw.bern');
child(2).b = readCoeffs2('bse.bern');
child(3).b = readCoeffs2('bnw.bern');
child(4).b = readCoeffs2('bne.bern');


box = importdata('obb.dat');
ctr = box(1,:)';
rng = box(2,:)';
axe = box(3:5,:)';





cl = colorcet( 'I2', 'N', 5 );


n = 200;
t = linspace(0,1,n)';

for i = 1:3
    figure;
    hold on
    
%     plot3( 0, 0, 0, 'k*' );
    
    switch i
        case 1
            surf_chebyshev( c, 1 );
        case 2
            s = eval_Bezier_patch_tpgrid( b, t, t );
            surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                'edgecolor','none','facecolor',cl(1,:),'diffusestrength',0.7, ...
                'facelighting','gouraud','specularstrength',0.15);
            
            h = plotOBB( ctr, rng, axe );
            set( h, 'facecolor', 'none', 'edgecolor', 'r' );
        case 3
            for j = 1:4
                s = eval_Bezier_patch_tpgrid( child(j).b, t, t );
                surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                    'edgecolor','none','facecolor',cl(1+j,:),'diffusestrength',0.7, ...
                    'facelighting','gouraud','specularstrength',0.15);
                boxb = OBB_Bernstein2( child(j).b );
                h = plotOBB( boxb.center, boxb.extents, boxb.axes );
                set( h, 'facecolor', 'none', 'edgecolor', 0.5*cl(1+j,:) );
            end
    end
    axis image vis3d
    view(3)
    camlight(30,30);
    camproj('persp');
end

