clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');


cl = colorcet( 'I2', 'N', 5 );
n = 200;
t = linspace(0,1,n)';
u = linspace(-1,1,n)';

% uvbox = [
%     -1.0000000000000000      -0.37073269459372837       -1.0000000000000000        9.8577654161285146E-003
%     -1.0000000000000000      -0.55555534117789551       -1.0000000000000000        1.9698044316964931E-002];
% uvbox = [-1.0000000000000000        1.0000000000000000       -1.0000000000000000        1.0000000000000000
%     -1.0000000000000000      -0.86602540378798298       -1.0000000000000000      -0.86602540378798298];
uvbox = [
    -1.0000000000000000        1.0000000000000000       -1.0000000000000000       -8.0000000001980986E-003
  -1.0000000000000000        1.0000000000000000       -1.0000000000000000       -1.9680784779652072E-013];

for isurf = 1:2
    
    
    
    for ipoly = 1:2
        if ipoly == 1
            root_c = readCoeffs2(sprintf('c_%d.cheb',isurf));
            region_b = readCoeffs2(sprintf('region_%d_x.bern',isurf));
            for ichild = 1:4
                child(ichild).b = readCoeffs2(sprintf('child_%d_%d_x.bern',isurf,ichild));
            end
        else
            root_c = readCoeffs2(sprintf('pn_%d.cheb',isurf));
            region_b = readCoeffs2(sprintf('region_%d_pn.bern',isurf));
            for ichild = 1:4
                child(ichild).b = readCoeffs2(sprintf('child_%d_%d_pn.bern',isurf,ichild));
            end
        end
        
        region_c = chgvar2( root_c, reshape(uvbox(isurf,:), 2, 2 ) );
        
        figure;
        hold on
        
        %         surf_chebyshev( root_c, 1 );
        
        [s, ~, ~, ~, lc ] = surf_chebyshev( region_c, 1, n, 0, 0, 1 );
        set( s, 'facecolor', 'none' );
        set( lc, 'color', 'g', 'linewidth', 1.5, 'linestyle', ':' );
        
        s = eval_Bezier_patch_tpgrid( region_b, t, t );
        lc = contourSurf(s(:,:,1),s(:,:,2),s(:,:,3));
        set( lc, 'linewidth', 1.5, 'linestyle', '--' );
        
        for ichild = 1:4
            s = eval_Bezier_patch_tpgrid( child(ichild).b, t, t );
            surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
                'edgecolor','none','facecolor',cl(ichild, :),'diffusestrength',0.7, ...
                'facelighting','gouraud','specularstrength',0.15);
        end
        
        axis image vis3d
        view(3)
        camlight(30,30);
        camproj('persp');
    end
    
end
