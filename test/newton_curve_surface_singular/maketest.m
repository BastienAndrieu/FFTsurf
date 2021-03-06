clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

EPS = 10.0 * eps('double');

format long

S.c = readCoeffs2( '/stck/bandrieu/Bureau/coeffstest/C1_test06.txt' );
S.c( abs(S.c) < EPS ) = 0;

S.d = cheb_diff2( S.c );
d2u = cheb_diff2( S.d(:,:,:,1) );
d2v = cheb_diff2( S.d(:,:,:,2) );
S.d2 = cat( 4, d2u, d2v(:,:,:,2) );

[Ms,Ns,~] = size( S.c );


W.c = [ 0.3, 0.2;
    0.0, 0.6;
    0.2, 0.0];

Nw = size( W.c, 1 );
t0 = 2*rand - 1
% t0 = 1 

Nd = 4;
D.c = ( 2*rand(Nd,3) - 1.0 ) ./ repmat( (1:Nd)'.^3, 1, 3 );

contact_order = 1;

Ng = max( (Ms + Ns) * Nw, Nd + contact_order + 2 );
f = 1;
gc = ICT2unstr( S.c, chebval1( W.c, cglpoints( 1, -1, Ng )' ) );
if contact_order >= 0
tc = cglpoints( 1, -1, Ng )';
gc = gc + ifcht( D.c, Ng ) .* repmat( f*( tc - t0 ).^(1+contact_order), 1, 3 );
end
G.c = fcht( gc );
G.d = cheb_diff1( G.c );
G.d2 = cheb_diff1( G.d );


% G.c

% G.c(1,:) = G.c(1,:) + 0.1 * (2*rand(1,3) - 1);






%%
writeCoeffs2( S.c, 'surf.cheb' );
writeCoeffs1( G.c, 'curv.cheb' );


%%
t = t0;
uv = chebval1( W.c, t );

fid = fopen( 'tuv.dat', 'w' );
fprintf( fid, '%16.15e\n', t );
fprintf( fid, '%16.15e %16.15e\n', uv(1), uv(2) );
fclose(fid);

g = chebval1( G.c, t )';
gt = chebval1( G.d, t )';
gtt = chebval1( G.d2, t )';

s = ICT2unstr( S.c, uv )';
su = ICT2unstr( S.d(:,:,:,1), uv )';
sv = ICT2unstr( S.d(:,:,:,2), uv )';
suu = ICT2unstr( S.d2(:,:,:,1), uv )';
suv = ICT2unstr( S.d2(:,:,:,2), uv )';
svv = ICT2unstr( S.d2(:,:,:,3), uv )';

r = s - g;
H = [
    gt'*gt - gtt'*r, -gt'*su        , -gt'*sv        ;
    -gt'*su        , su'*su + suu'*r, su'*sv + suv'*r;
    -gt'*sv        , su'*sv + suv'*r, sv'*sv + svv'*r]
rank(H)

J = [-gt, su, sv]
rank(J)

% return


tuv = [t,uv]

g = chebval1( G.c, t );
s = ICT2unstr( S.c, uv );
xyz = 0.5 * ( g + s )
% return

% g = chebval1( G.c, t )';
g1 = chebval1( G.d, t )';
g2 = chebval1( G.d2, t )';

s = ICT2unstr( S.c, uv )';
s1 = zeros(3,2);
s2 = zeros(3,3);
for ivar = 1:2
    s1(:,ivar) = ICT2unstr( S.d(:,:,:,ivar), uv )';
end
for ivar = 1:3
    s2(:,ivar) = ICT2unstr( S.d2(:,:,:,ivar), uv )';
end

n = cross( s1(:,1), s1(:,2) );

wt = s1 \ g1;

y = g2 - ( wt(1)^2 * s2(:,1) + 2 * wt(1)*wt(2) * s2(:,2) + wt(2)^2 * s2(:,3) );



fprintf( '|gt.n| = %e\n |y.n| = %e\n', abs(dot(g1,n)), abs(dot(y,n)) );





%% plots
n = 100;
t = linspace(-1,1,n)';


figure;
hold on

surf_chebyshev( S.c, 1 );

g = chebval1( G.c, t );
g0 = chebval1( G.c, t0 );

plot3( g(:,1), g(:,2), g(:,3), 'k-' );
plot3( g0(1), g0(2), g0(3), 'k*' );

axis image vis3d
view(3)
camproj( 'persp' );
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );