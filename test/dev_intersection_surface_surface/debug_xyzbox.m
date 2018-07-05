clc; clear; close all



cl = colorcet( 'I2', 'N', 2 );
cls = CC( cl, 0.0, 0.8, 1.5 );
clb = CC( cl, 0.0, 1.5, 0.5 );

figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cls(isurf,:), 'specularstrength', 0 );
    
    box = importdata( sprintf('xyzbox%d.dat',isurf) );
    ctr = box(1,:)';
    rng = box(2,:)';
    axe = box(3:5,:)';
    
    h = plotOBB( ctr, rng, axe );
    set( h, 'facecolor', 'none', 'edgecolor', clb(isurf,:), 'linestyle', '--' );
end


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );