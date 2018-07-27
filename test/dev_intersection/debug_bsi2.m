clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

cl = colorcet( 'I2', 'N', 2 );

bc = readCoeffs1('debugbsi_c.bern');
bs = readCoeffs2('debugbsi_s.bern');

cc = readCoeffs1('debugbsi_c.cheb');
cs = readCoeffs2('debugbsi_s.cheb');

tuv = [
    0.748620759138596 -0.178813934326172e-6  -0.178813934326172e-6
%     -0.83333333333333326      -0.93333333333333335      -0.91666666666666663
    ];

n = 200;
t = linspace( -1, 1, n )';

figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

s = surf_chebyshev(cs, 1, n);
set(s, 'specularstrength', 0);

c = chebval1(cc, t);
plot3( c(:,1), c(:,2), c(:,3), 'k-' );

% s = eval_Bezier_patch_tpgrid( bs, t, t );
% surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
%     'edgecolor','none','facecolor','w','diffusestrength',0.7, ...
%     'facelighting','gouraud','specularstrength',0.15);
%
% c = evalBernstein( bc, t );
% plot3( c(:,1), c(:,2), c(:,3), 'k-' );

s = icht2_pc(cs, tuv(:,2:3));
plot3( s(:,1), s(:,2), s(:,3), 'ro' );

c = chebval1(cc, tuv(:,1));
plot3( c(:,1), c(:,2), c(:,3), 'b.' );

axis image vis3d
view(3)
camproj('persp');
camlight(30,30);

%%
ds = cheb_diff2(cs);
d2u = cheb_diff2(ds(:,:,:,1));
d2v = cheb_diff2(ds(:,:,:,2));
d2s = cat(4, d2u, d2v(:,:,:,2));

dc = cheb_diff1(cc);
d2c = cheb_diff1(dc);

g1 = chebval1( dc, tuv(1) )';
g2 = chebval1( d2c, tuv(1) )';

s1 = zeros(3,2);
s2 = zeros(3,3);
for ivar = 1:2
    s1(:,ivar) = ICT2unstr( ds(:,:,:,ivar), tuv(2:3) )';
end
for ivar = 1:3
    s2(:,ivar) = ICT2unstr( d2s(:,:,:,ivar), tuv(2:3) )';
end

n = cross( s1(:,1), s1(:,2) );
    n = n / norm(n);
    
    wt = s1 \ g1;
    
    y = g2 - ( wt(1)^2 * s2(:,1) + 2 * wt(1)*wt(2) * s2(:,2) + wt(2)^2 * s2(:,3) );
    
    fprintf( '|gt.n|/|gt| = %e\n |y.n|/|y| = %e\n', abs(dot(g1,n))/norm(g1), abs(dot(y,n))/norm(y) );
    fprintf( '     |gt.n| = %e\n     |y.n| = %e\n', abs(dot(g1,n))         , abs(dot(y,n))         );