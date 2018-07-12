clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');

cl = colorcet( 'I2', 'N', 2 );
cl = CC( cl, 0.0, 0.8, 1.5 );

suffix = {'x';'pn'};

for i = 1:2
    
    figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
    hold on
    
    for isurf = 1:2
        c = readCoeffs2( sprintf('surfroot%d_%s.cheb', isurf, suffix{i} ) );
        b = readCoeffs2( sprintf('surfroot%d_%s.bern', isurf, suffix{i} ) );
        
        si = surf_chebyshev( c, 1 );
        set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
        
        surf( b(:,:,1), b(:,:,2), b(:,:,3), 'facecolor', 'none', 'edgecolor', 0.5*cl(isurf,:) );
    end
    
    
    axis image vis3d
    view(3)
    camproj('persp');
    % camlight(30,30);
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
    
    
end
