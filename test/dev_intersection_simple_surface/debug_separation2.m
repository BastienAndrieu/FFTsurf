clc; clear; close all


addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');




ts = 0.510797450327141;
us = 0.751788519744514;
vs = 0.776852809891227;


bc = readCoeffs1('root_c_bezier.bern');
bs = readCoeffs2('root_s_bezier.bern');



[ bc, ~ ] = subdivide_deCasteljau( bc, ts );
[ ~, bs, ~, ~  ] = subdivide_Bezier_patch( bs, [us,vs] );



xyz = bc(end,:);



for i = 1:3
    bc(:,i) = bc(:,i) - xyz(i);
    bs(:,:,i) = bs(:,:,i) - xyz(i);
end

xc = bc;
xs = reshape( bs, [], 3 );

L = sqrt( max( [ max(sum(xc.^2,2)), max(sum(xs.^2,2)) ] ) );

EPSsqr = (1e-9)^2;
ic = find( sum( xc.^2, 2 ) > EPSsqr );
is = find( sum( xs.^2, 2 ) > EPSsqr );

yc = xc(ic,:);
ys = xs(is,:);


yc = yc ./ repmat( sqrt( sum( yc.^2, 2 ) ), 1, 3 )
ys = ys ./ repmat( sqrt( sum( ys.^2, 2 ) ), 1, 3 )


R = [
    0.71294066486393737       -1.2369239509165280E-002 -0.70100970004788243     
   1.7349605821020062E-002  0.99969899117785521        0.0000000000000000     
  0.70100970004788243       -1.2162241972542268E-002  0.71294066486393737];
yc = yc * R;
ys = ys * R;

[vec,stat] = separating_plane( yc, ys );

vec = R*vec;

bc = bc / L;
bs = bs / L;



% figure;
% hold all
% for i = 1:3
%     quiver3( 0,0,0,R(i,1),R(i,2),R(i,3),0);
% end
% axis image
% view(3)




n = 300;
t = linspace( 0, 1, n )';


figure;
hold on

plot3( 0, 0, 0, 'k*', 'markersize', 12 );

s = eval_Bezier_patch_tpgrid( bs, t, t );
surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
    'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
    'facelighting','gouraud','specularstrength',0.15);
surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), ...
    'facecolor', 'none', 'edgecolor', 'r', 'linestyle', '-' );


c = evalBernstein( bc, t );
plot3( c(:,1), c(:,2), c(:,3), 'k-' );
plot3( bc(:,1), bc(:,2), bc(:,3), 'b.-' );

if stat == 0
    [ q, p ] = plot_separating_plane( vec );
end

plot3( yc(:,1), yc(:,2), yc(:,3), 'b*' );
plot3( ys(:,1), ys(:,2), ys(:,3), 'r*' );

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);