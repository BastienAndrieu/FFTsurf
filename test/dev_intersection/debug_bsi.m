clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');


cl = colorcet( 'I2', 'N', 2 );

uvbox = [
    0.0000000000000000        1.0000000000000000       -1.0000000000000000        0.0000000000000000     
  -1.0000000000000000        0.0000000000000000        0.0000000000000000        1.0000000000000000
%     -1.0000000000000000        1.0000000000000000       -1.0000000000000000        1.0000000000000000     
%   -1.0000000000000000        1.0000000000000000       -1.0000000000000000        1.0000000000000000
    ];

tbox = [0,1];
ubox = [-1,-1;0,0];

[icurv, ivar, ival] = deal(1,2,1);
jvar = 1 + mod(ivar,2);

tuvxyz = [
    -1.0000000000000000       -1.0000000000000000        1.0000000000000000 1.0000000000000000      -0.11755705045849460       0.16180339887498874
];

uv = zeros(2,2,size(tuvxyz,1));
for i = 1:size(tuvxyz,1)
    uv(ivar,icurv) = uvbox(icurv,2*(ivar-1)+ival);
    uv(jvar,icurv) = uvbox(icurv,2*jvar-1) + 0.5 * ( uvbox(icurv,2*jvar) - uvbox(icurv,2*jvar-1) ) * (tuvxyz(i,1)+1);
    uv(:,1+mod(icurv,2)) = tuvxyz(i,2:3)';
end
uv

% uv = reshape( [0.0000000000000000      -0.69176322065777784      -0.69176322065777784      -0.69176322065777784], [2,2,1] );


uvxyz = [
    4.3481470779374379E-002  0.50000000000000000       0.54348147077937448       0.45651852922062569       -3.7777356481573335E-002 -0.43440752813142158       0.62504137348464284
    4.3481469154357910E-002  0.50000000000000000       0.54348146915435791       0.45651853084564209       -3.7777355069734557E-002 -0.43440752813142158       0.62504137255360726
    1.3080661123996255E-002  0.51308066112399620       0.52616132224799261       0.50000000000000000       -1.1364675330399990E-002 -0.44577220346182156       0.60778005247054423
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
    cr = readCoeffs2( sprintf('surfroot%d_x.cheb', isurf) );
    
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
        si = surf_chebyshev( c, 1 );
        set( si, 'facecolor', cl(isurf,:), 'specularstrength', 0 );
    end
    
    %     b = readCoeffs2( sprintf('debugbsi_reg%d.bern', isurf) );
    %     surf( b(:,:,1), b(:,:,2), b(:,:,3), 'facecolor', 'none', 'edgecolor', 0.25*cl(isurf,:) );
    
    %     xi = ICT2unstr(cr, uvxyz(:,2*isurf-1:2*isurf));
    %     plot3( xi(:,1), xi(:,2), xi(:,3), 'x', 'color', 0.4*cl(isurf,:) );
    
    %     uvi = reshape( uv(:,isurf,:), 2,[] )';
    %     yi = ICT2unstr(cr, uvi);
    %     plot3( yi(:,1), yi(:,2), yi(:,3), '+', 'color', 0.4*cl(isurf,:) );
    
    if isurf == icurv
        c = border_parameterization( c, ivar, ival );
        G.c = c;
        %         c = chgvar1(c, tbox);
        %         writeCoeffs1(c,'../newton_curve_surface_singular/curv.cheb');
        
        g = chebval1( c, linspace(-1,1,100)' );
        plot3( g(:,1), g(:,2), g(:,3), '-', 'color', 0.4*cl(isurf,:) );
        
        z = chebval1(c, tuvxyz(:,1));
    else
        z = ICT2unstr(cr, tuvxyz(:,2:3));
    end
    
    plot3( z(:,1), z(:,2), z(:,3), 'o', 'color', 0.4*cl(isurf,:) );
end





% plot3( uvxyz(:,5), uvxyz(:,6), uvxyz(:,7), 'k.', 'markersize', 8 );
% plot3( tuvxyz(:,4), tuvxyz(:,5), tuvxyz(:,6), 'k.', 'markersize', 8 );


axis image vis3d
view(3)
camproj('persp');
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );


%%
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
EPS = 1e-9;
cr = chgvar2( S.c, reshape(uvbox(1+mod(icurv,2),:),2,2) );
S.b = chebyshev2bezier_2(cr);
G.b = chebyshev2bezier_1(G.c);

seps = reshape(S.b, [], 3);
seps = seps - repmat(tuvxyz(1,4:6),size(seps,1),1);
seps(sqrt(sum(seps.^2,2)) < EPS,:) = [];
seps = seps + repmat(tuvxyz(1,4:6),size(seps,1),1);


sepc = G.b;
sepc = sepc - repmat(tuvxyz(1,4:6),size(sepc,1),1);
sepc(sqrt(sum(sepc.^2,2)) < EPS,:) = [];
sepc = sepc + repmat(tuvxyz(1,4:6),size(sepc,1),1);


[ vec, stat ] = separating_plane( sepc, seps )




qs = myconvhulln(seps);
qc = myconvhulln(sepc);

figure;
hold on

for isurf = 1:2
    if isurf == icurv
        q = qc;
        sep = sepc;
    else
        q = qs;
        sep = seps;
    end
    trisurf(q, sep(:,1), sep(:,2), sep(:,3), ...
        'facecolor', cl(isurf,:), 'edgecolor', 'none' );
end

if stat == 0
    plot_separating_plane(vec);
end

axis image vis3d
view(120,30)
camproj('persp');
camlight(30,30);


