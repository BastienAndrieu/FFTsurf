clc; clear; close all


addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
% 
% tbox = [-0.87424144808651783      -0.49166357681374367];
% uvbox = [-0.94999999999999996      -0.90000000000000002      -0.78681730327786603      -0.49249721913236932];

% bc = readCoeffs1('region_c_bezier.bern');
% bs = readCoeffs2('region_s_bezier.bern');

% A = C2Bmatrix( size(bc,1) );
% cc = A \ bc;
% cc = chgvar1( cc, tbox );
% bc = A * cc;
% 
% cs = readCoeffs2('c_1.cheb');
% cs = chgvar2( cs, reshape(uvbox,2,2) );
% bs = chebyshev2bezier_2( cs );

bc = readCoeffs1( 'region_c_bezier.bern' );
bs = readCoeffs2( 'region_s_bezier.bern' );
% xyz = [1.0000000000000000      -0.63111829381363049       0.12277330464232007];
% xyz = [-0.90000000000000002        1.5005561737555561E-002   0.0000000000000000];
xyz = importdata('debug_separation_xyz.dat');


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
% format long
% sqrt( sum( yc.^2, 2 ) )
% sqrt( sum( ys.^2, 2 ) )
% format short

yc = yc ./ repmat( sqrt( sum( yc.^2, 2 ) ), 1, 3 )
ys = ys ./ repmat( sqrt( sum( ys.^2, 2 ) ), 1, 3 )

% R = eye(3);
% w = rand(3,1);
% w = w / norm(w);
% [u,v] = completeBase(w);
% R = [u,v,w];
% R = [
%     -0.084837363255110   0.929144767855911   0.359850833209964
%   -0.213207351598671  -0.369716107798923   0.904352046969705
%    0.973316622184272                   0   0.229466234944924];
% 
R = [
    0.37552734040600005       -7.2700052163319775E-002 -0.90636657096156981     
  0.19359456513797468       0.96252114434903846       -1.3877787807814457E-017
  0.90636657096156992      -0.17546764216090241       0.37552734040600005];



yc = yc * R;
ys = ys * R;




[vec,stat] = separating_plane( yc, ys );
stat
vec = R*vec;
yc = yc * R';
ys = ys * R';

n = 100;
t = linspace( 0, 1, n )';

bc = bc / L;
bs = bs / L;


figure;
hold on

plot3( 0, 0, 0, 'k*' );

s = eval_Bezier_patch_tpgrid( bs, t, t );
surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
    'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
    'facelighting','gouraud','specularstrength',0.15);
surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), 'facecolor', 'none', 'edgecolor', 'r' );

c = evalBernstein( bc, t );
plot3( c(:,1), c(:,2), c(:,3), 'k-' );
plot3( bc(:,1), bc(:,2), bc(:,3), 'b.-' );

if stat == 0
    [ q, p ] = plot_separating_plane( vec );
end

plot3( yc(:,1), yc(:,2), yc(:,3), 'b*' );
plot3( ys(:,1), ys(:,2), ys(:,3), 'r*' );

axis image
view(3)
camlight(30,30);


return
figure;
hold on


plot3( 0, 0, 0, 'k*' );

plot3( yc(:,1), yc(:,2), yc(:,3), 'b*' );
plot3( ys(:,1), ys(:,2), ys(:,3), 'r*' );

if stat == 0
    [ q, p ] = plot_separating_plane( vec );
end

axis image
view(3)



return
%%

n = 100;
t = linspace(0,1,n)';


figure;
hold on

s = eval_Bezier_patch_tpgrid( bs, t, t );
surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
    'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
    'facelighting','gouraud','specularstrength',0.15);
% surf( bs(:,:,1), bs(:,:,2), bs(:,:,3), ...
%     'facecolor', 'none', 'edgecolor', 'r', 'linestyle', '--' );


c = evalBernstein( bc, t );
plot3( c(:,1), c(:,2), c(:,3), 'k-' );
% plot3( bc(:,1), bc(:,2), bc(:,3), 'b.--' );

% plot3( xyz(1), xyz(2), xyz(3), 'k*', 'markersize', 4 );

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);


