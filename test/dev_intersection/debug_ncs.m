clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');


n = 200;
t = linspace(0,1,n)';

bc = readCoeffs1('debugncs_regc.bern');
bs = readCoeffs2('debugncs_regs.bern');


cc = readCoeffs1('debugncs_rootc.cheb');
cs = readCoeffs2('debugncs_roots.cheb');

tbox = importdata('debugncs_uvboxc.dat');
cc = chgvar1(cc, tbox);
uvbox = importdata('debugncs_uvboxs.dat');
cs = chgvar2(cs, reshape(uvbox, 2, 2));

writeCoeffs1(cc,'../newton_curve_surface_singular/curv.cheb');
writeCoeffs2(cs,'../newton_curve_surface_singular/surf.cheb');

%%
if 0
    figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
    hold on
    
    s = eval_Bezier_patch_tpgrid(bs, t, t);
    surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
        'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
        'facelighting','gouraud','specularstrength',0.15);
    
    c = evalBernstein(bc, t);
    plot3( c(:,1), c(:,2), c(:,3), 'k-' );
    
    s0 = eval_Bezier_patch_tpgrid(bs, 0.5, 0.5);
    plot3( s0(1), s0(2), s0(3), 'ro' );
    
    c0 = evalBernstein(bc, 0.5);
    plot3( c0(1), c0(2), c0(3), 'b*' );
    
    axis image vis3d
    view(3)
    camproj('persp');
    % camlight(30,30);
    [azl,ell] = deal( 120, 60 );
    [xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
    light( 'style', 'infinite', 'position', [xl,yl,zl] );
    light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
end


%%
figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

surf_chebyshev( cs, 1 );

g = chebval1( cc, linspace(-1,1,n)' );
plot3( g(:,1), g(:,2), g(:,3), 'k-' );

surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), 'facecolor', 'none', 'edgecolor', 'r' );
plot3( bc(:,1), bc(:,2), bc(:,3), 'b.-');

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);