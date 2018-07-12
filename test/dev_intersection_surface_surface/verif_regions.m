clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');

uvbox = [
    -1.0000000000000000        1.0000000000000000       -1.0000000000000000       -8.0000000001874683E-003
    -1.0000000000000000        1.0000000000000000       -1.0000000000000000       -1.9061835443423547E-013
    ];
cl = colorcet( 'I2', 'N', 2 );
% cl = linspecerp(2);
cl = CC( cl, 0.0, 0.8, 1.5 );


figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

for isurf = 1:2
    c = readCoeffs2( sprintf('c_%d.cheb', isurf) );
    c = chgvar2( c, reshape(uvbox(isurf,:),2,2) );
    
    
    si = surf_chebyshev( c, 1 );
    set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    
    %         b = chebyshev2bezier_2( c );
    %     quiver3( b(1,1,1), b(1,1,2), b(1,1,3), ...
    %         b(2,1,1)- b(1,1,1), b(2,1,2)- b(1,1,2), b(2,1,3)- b(1,1,3), ...
    %         'r' );
    %     quiver3( b(1,1,1), b(1,1,2), b(1,1,3), ...
    %         b(1,2,1)- b(1,1,1), b(1,2,2)- b(1,1,2), b(1,2,3)- b(1,1,3), ...
    %         'g' );
end



icurv = 1;
ivar = 1;
ival = 1;

c = readCoeffs2( sprintf('c_%d.cheb', icurv) );
c = chgvar2( c, reshape(uvbox(icurv,:),2,2) );
c = border_parameterization( c, ivar, ival );

n = 200;
u = linspace(-1,1,n)';
g = chebval1( c, u );
plot3( g(:,1), g(:,2), g(:,3), 'k-' );



axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


%%
tbox = [0.0000000000000000        1.0000000000000000];
uvbox = [0.0000000000000000        1.0000000000000000      -0.50000000000009526       -1.9061835443423547E-013];
isurf = 1 + mod(icurv,2);


cc = chgvar1( c, tbox );
cs = readCoeffs2( sprintf('c_%d.cheb', isurf) );
cs = chgvar2( cs, reshape(uvbox,2,2) );

bc = chebyshev2bezier_1( cc );
bs = chebyshev2bezier_2( cs );



figure( 'units', 'normalized', 'position',[.15 .15 .7 .7 ] );
hold on

si = surf_chebyshev( cs, 1 );
set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );

surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), 'facecolor', 'none', 'edgecolor', 0.4*cl(isurf,:) );

g = chebval1( cc, u );
plot3( g(:,1), g(:,2), g(:,3), 'k-' );

plot3( bc(:,1), bc(:,2), bc(:,3), '.-', 'color', 0.4*cl(icurv,:) );

axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


%%
x = bc(end,:);

sepc = bc - repmat( x, size(bc,1), 1 );
seps = reshape(bs,[],3) - repmat( x, size(bs,1)*size(bs,2), 1 );

tol = 1e-9;
sepc( sqrt(sum(sepc.^2,2)) < tol, : ) = [];
seps( sqrt(sum(seps.^2,2)) < tol, : ) = [];
sepc = sepc ./ repmat( sqrt(sum(sepc.^2,2)), 1, 3 );
seps = seps ./ repmat( sqrt(sum(seps.^2,2)), 1, 3 );

[xsp, ysp, zsp] = sphere(100);

figure;
hold on

plot3( sepc(:,1), sepc(:,2), sepc(:,3), '.', 'color', 'b' );
plot3( seps(:,1), seps(:,2), seps(:,3), '.', 'color', 'r' );

surf( xsp, ysp, zsp, 'facecolor', 'w', 'edgecolor', 'none', 'facelighting',' phong', 'facealpha', 0.2 );

axis image vis3d
view(3)
camlight(30,30);

addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
[vec,stat] = separating_plane( sepc, seps )


