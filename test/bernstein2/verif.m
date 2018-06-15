clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');


%%
c = readCoeffs2('c.cheb');

b = readCoeffs2('b.bern');

child(1).b = readCoeffs2('bsw.bern');
child(2).b = readCoeffs2('bse.bern');
child(3).b = readCoeffs2('bnw.bern');
child(4).b = readCoeffs2('bne.bern');

child_only_u(1).b = readCoeffs2('bw.bern');
child_only_u(2).b = readCoeffs2('be.bern');

child_only_v(1).b = readCoeffs2('bs.bern');
child_only_v(2).b = readCoeffs2('bn.bern');


%%
n = 200;
t = linspace(0,1,n)';

cl = colorcet( 'I2', 'N', 9 );
cl = cl(randperm(size(cl,1)),:);

for i = 1:5
    figure;
    hold on
    
    switch i
        case 1
            surf_chebyshev( c, 1, n, 11 );
            
        case 2
            s = eval_Bezier_patch_tpgrid( b, t, t );
            surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                'edgecolor','none','facecolor',cl(1,:),'diffusestrength',0.7, ...
                'facelighting','gouraud','specularstrength',0.15);
            
        case 3
            for j = 1:2
                s = eval_Bezier_patch_tpgrid( child_only_u(j).b, t, t );
                surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                    'edgecolor','none','facecolor',cl(1+j,:),'diffusestrength',0.7, ...
                    'facelighting','gouraud','specularstrength',0.15);
            end
            
        case 4
            for j = 1:2
                s = eval_Bezier_patch_tpgrid( child_only_v(j).b, t, t );
                surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                    'edgecolor','none','facecolor',cl(1+2+j,:),'diffusestrength',0.7, ...
                    'facelighting','gouraud','specularstrength',0.15);
            end
            
        case 5
            for j = 1:4
                s = eval_Bezier_patch_tpgrid( child(j).b, t, t );
                surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                    'edgecolor','none','facecolor',cl(1+2+2+j,:),'diffusestrength',0.7, ...
                    'facelighting','gouraud','specularstrength',0.15);
            end
            
    end
   
    axis image vis3d
    view(3)
    camlight(30,30);
    camproj('persp');
end