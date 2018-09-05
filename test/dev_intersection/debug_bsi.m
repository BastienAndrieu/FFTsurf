clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');

cl = colorcet( 'I2', 'N', 2 );

uvbox = [
    -1.0000000000000000        1.0000000000000000       -1.0000000000000000        1.0000000000000000     
  -1.0000000000000000      -0.98698748125962998       -1.0000000000000000       -6.3451034861136291E-002
% 0.0000000000000000        1.5258789062500000E-005 -0.30496215820312500      -0.30494689941406250     
%  -0.10777282714843750      -0.10775756835937500       -4.9316406250000000E-002  -4.9301147460937500E-002
%     -3.0575494019804609E-005  0.49997244140091862       0.49994384986693124        1.0000000000000000     
%   -1.0000000000000000       -2.7497621712598644E-005  -1.0000000000000000       -2.7457768678921469E-006
    ];

tbox = [
    -1.0000000000000000        1.0000000000000000
    ];
ubox = [
    -1.0000000000000000      -0.98698748125962998       -1.0000000000000000       -6.3451034861136291E-002
%     -1.0000000000000000        0.0000000000000000       -1.0000000000000000        0.0000000000000000
%     -0.10777282714843750      -0.10775756835937500       -4.9316406250000000E-002  -4.9301147460937500E-002
% -3.0575494044112767E-005  -2.7497621712598644E-005  -1.0000000000000000       -2.7457768678921469E-006
];

[icurv, ivar, ival] = deal(1,2,2);
jvar = 1 + mod(ivar,2);

tuvxyz = [
    0.39010773412313421      -0.10776183717734601       -4.9314229082246912E-002  -4.7846698502724474E-017   1.4877639051280589E-014   1.6028844918025698E-015
%     -1.0000000000000000      -0.50002376865584797       -1.0000000000000224       0.67744753604335384       0.54810455295093197       0.56146345939272302
%     -1.0000000000000000       -3.0575494044112767E-005 -0.99999999999999856       0.68401705066906149       0.52844689545221724       0.53512981220642375 
];

uv = zeros(2,2,size(tuvxyz,1));
for i = 1:size(tuvxyz,1)
    uv(ivar,icurv,i) = uvbox(icurv,2*(ivar-1)+ival);
    uv(jvar,icurv,i) = uvbox(icurv,2*jvar-1) + 0.5 * ( uvbox(icurv,2*jvar) - uvbox(icurv,2*jvar-1) ) * (tuvxyz(i,1)+1);
    uv(:,1+mod(icurv,2),i) = tuvxyz(i,2:3)';
end
uv

% uv = reshape( [0.0000000000000000      -0.69176322065777784      -0.69176322065777784      -0.69176322065777784], [2,2,1] );


uvxyz = [
    -0.50002376865578813       0.99999999999996714      -0.50002376865584797       -1.0000000000000224 0.67744753604335384       0.54810455295093197       0.56146345939272302
%     -3.0575494019804609E-005  0.99999999999998512       -3.0575494044112767E-005 -0.99999999999999856 0.68401705066906149       0.52844689545221724       0.53512981220642375
    ];

% for i = 1:size(tuvxyz,1)
%     reshape( uvxyz(i,1:4), 2, 2 )
% end

% t = uvxyz(2*(icurv-1)+jvar);
% t = -1 + 2 * (t - uvbox(icurv,2*jvar-1)) / (uvbox(icurv,2*jvar) - uvbox(icurv,2*jvar-1))
% erruv = uvxyz(2*(icurv-1)+ivar) - uvbox(icurv,2*(ivar-1)+ival)


figure('units', 'normalized', 'position',[.15 .15 .7 .7 ]);
hold on

for isurf = 1:2
    cr = readCoeffs2( sprintf('surfroot%2.2d_x.cheb', isurf) );
    
    if isurf ~= icurv
        S.c = cr;
    end
    
    if 0%isurf ~= icurv
        c = chgvar2( cr, ubox );
        writeCoeffs2(c,'../newton_curve_surface_singular/surf.cheb');
    else
        c = chgvar2( cr, reshape(uvbox(isurf,:),2,2) );
    end
    
    if 1%isurf ~= icurv
%         c = chgvar2( cr, reshape(ubox,2,2) );
        si = surf_chebyshev( c, 1 );
        set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
        
%         x = cat(3, ifcht2d(c(:,:,1)), ifcht2d(c(:,:,2)),  ifcht2d(c(:,:,3)));
%         plot3(x(:,:,1), x(:,:,2), x(:,:,3), '.', 'color', 0.6*cl(isurf,:));
    end
    
    %     b = readCoeffs2( sprintf('debugbsi_reg%d.bern', isurf) );
    %     surf( b(:,:,1), b(:,:,2), b(:,:,3), 'facecolor', 'none', 'edgecolor', 0.25*cl(isurf,:) );
    
    %     xi = ICT2unstr(cr, uvxyz(:,2*isurf-1:2*isurf));
    %     plot3( xi(:,1), xi(:,2), xi(:,3), 'x', 'color', 0.4*cl(isurf,:) );
    
%         uvi = reshape( uv(:,isurf,:), 2,[] )';
%         yi = ICT2unstr(cr, uvi);
%         plot3( yi(:,1), yi(:,2), yi(:,3), '+', 'color', 0.4*cl(isurf,:) );
    
    if isurf == icurv
        c = border_parameterization( c, ivar, ival );
        G.c = c;
                c = chgvar1(c, tbox);
        %         writeCoeffs1(c,'../newton_curve_surface_singular/curv.cheb');
        
        g = chebval1( c, linspace(-1,1,100)' );
        plot3( g(:,1), g(:,2), g(:,3), '-', 'color', 0.2*cl(isurf,:), 'linewidth', 1.5 );
        
        z = chebval1(c, tuvxyz(:,1));
    else
        z = ICT2unstr(cr, tuvxyz(:,2:3));
    end
    
%     plot3( z(:,1), z(:,2), z(:,3), 'o', 'color', 0.4*cl(isurf,:) );
end





% plot3( uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), 'k.', 'markersize', 8 );
plot3( tuvxyz(:,4), tuvxyz(:,5), tuvxyz(:,6), 'k.', 'markersize', 8 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );



%%
figure;
for isurf = 1:2
    subplot(1,2,isurf);
    hold on
    rectangle('position',[uvbox(isurf,[1,3]), uvbox(isurf,[2,4]) - uvbox(isurf,[1,3])]);
    plot(uv(1,isurf,:), uv(2,isurf,:), '*');
    axis image
end


%%
figure;
hold on

cr = chgvar2(S.c, reshape(ubox,2,2));
bs = chebyshev2bezier_2(cr);

isurf = 1 + mod(icurv,2);
si = surf_chebyshev( cr, 1 );
set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );


cr = chgvar1(G.c, tbox);
bc = chebyshev2bezier_1(cr);

g = chebval1( cr, linspace(-1,1,100)' );
plot3( g(:,1), g(:,2), g(:,3), '-', 'color', 0.2*cl(icurv,:), 'linewidth', 1 );

boxc = OBB_Bernstein1(bc);
boxs = OBB_Bernstein2(bs);

b = plotOBB(boxc.center, boxc.extents, boxc.axes);
set(b, 'edgecolor', 0.5*cl(icurv,:), 'facecolor', 'none');

b = plotOBB(boxs.center, boxs.extents, boxs.axes);
set(b, 'edgecolor', 0.5*cl(isurf,:), 'facecolor', 'none');


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );





%%
return
S.d = cheb_diff2( S.c );
d2u = cheb_diff2( S.d(:,:,:,1) );
d2v = cheb_diff2( S.d(:,:,:,2) );
S.d2 = cat( 4, d2u, d2v(:,:,:,2) );

G.d = cheb_diff1( G.c );
G.d2 = cheb_diff1( G.d );

t = tuvxyz(:,1);
uv = tuvxyz(:,2:3);

for i = 1:size(t,1)
    g1 = chebval1( G.d, t(i) )';
    g2 = chebval1( G.d2, t(i) )';
    
    s1 = zeros(3,2);
    s2 = zeros(3,3);
    for ivar = 1:2
        s1(:,ivar) = ICT2unstr( S.d(:,:,:,ivar), uv(i,:) )';
    end
    for ivar = 1:3
        s2(:,ivar) = ICT2unstr( S.d2(:,:,:,ivar), uv(i,:) )';
    end
    
    n = cross( s1(:,1), s1(:,2) );
    n = n / norm(n);
    
    wt = s1 \ g1;
    
    y = g2 - ( wt(1)^2 * s2(:,1) + 2 * wt(1)*wt(2) * s2(:,2) + wt(2)^2 * s2(:,3) );
    
    fprintf( '|gt.n|/|gt| = %e\n |y.n|/|y| = %e\n', abs(dot(g1,n))/norm(g1), abs(dot(y,n))/norm(y) );
    fprintf( '     |gt.n| = %e\n     |y.n| = %e\n', abs(dot(g1,n))         , abs(dot(y,n))         );
    
    s = ICT2unstr( S.c, uv(i,:) )';
    g = chebval1( G.c, t(i) )';
    x = 0.5*(s + g);
    
    g1 = g1 / norm(g1);
    quiver3( x(1), x(2), x(3), n(1), n(2), n(3), 0.2, 'r' );
    quiver3( x(1), x(2), x(3), g1(1), g1(2), g1(3), 0.2, 'b' );
end


%%
% return
EPS = 1e-9;
% cr = chgvar2( S.c, reshape(uvbox(1+mod(icurv,2),:),2,2) );
cr = chgvar2( S.c, reshape(ubox,2,2) );
S.b = chebyshev2bezier_2(cr);
G.b = chebyshev2bezier_1(G.c);
% G.b = chebyshev2bezier_1( chgvar1(G.c, tbox) );

seps = reshape(S.b, [], 3);
seps = seps - repmat(tuvxyz(1,4:6),size(seps,1),1);
sqrt(sum(seps.^2,2))
seps(sqrt(sum(seps.^2,2)) < EPS,:) = [];
% seps = seps + repmat(tuvxyz(1,4:6),size(seps,1),1);


sepc = G.b;
sepc = sepc - repmat(tuvxyz(1,4:6),size(sepc,1),1);
sepc(sqrt(sum(sepc.^2,2)) < EPS,:) = [];
% sepc = sepc + repmat(tuvxyz(1,4:6),size(sepc,1),1);


seps = seps ./ repmat(sqrt(sum(seps.^2,2)), 1, 3);
sepc = sepc ./ repmat(sqrt(sum(sepc.^2,2)), 1, 3);

[ vec, stat ] = separating_plane( sepc, seps )




qs = myconvhulln(seps);
qc = myconvhulln(sepc);

figure;
hold on

plot3(0,0,0,'kx');

% surf(S.b(:,:,1)-tuvxyz(1,4), S.b(:,:,2)-tuvxyz(1,5), S.b(:,:,3)-tuvxyz(1,6), ...
%     'facecolor', 'none', 'edgecolor', 0.5*cl(1+mod(icurv,2),:));
% plot3(G.b(:,1)-tuvxyz(1,4), G.b(:,2)-tuvxyz(1,5), G.b(:,3)-tuvxyz(1,6), ...
%     '.-', 'color', 0.5*cl(icurv,:));

for isurf = 1:2
    if isurf == icurv
        q = qc;
        sep = sepc;
    else
        q = qs;
        sep = seps;
    end
%     plot3( sep(:,1), sep(:,2), sep(:,3), ...
%         '.', 'color', 0.7*cl(isurf,:));
%     trisurf(q, sep(:,1), sep(:,2), sep(:,3), ...
%         'facecolor', cl(isurf,:), 'edgecolor', 'none', 'facealpha', 1 );
[ fa, ed ] = patch_ridges( q, sep );
set(fa, 'facecolor', cl(isurf,:), 'specularstrength', 0.1);
end

if stat == 0
    plot_separating_plane(vec);
end

% axis image vis3d
axis( repmat([-1,1],1,3) );
daspect([1,1,1])
axis vis3d
view(120,30)
camproj('persp');
camlight(30,30);


