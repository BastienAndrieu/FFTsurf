clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Newton/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');

%% curve and surface
C.x = readCoeffs1('../dev_intersection/debugbsi_c.cheb');
% C.x = readCoeffs1( 'curv.cheb' );
C.xt = cheb_diff1( C.x );
C.xtt = cheb_diff1( C.x );

% S.x = readCoeffs2( 'surf.cheb' );
S.x = readCoeffs2('../dev_intersection/debugbsi_s.cheb');
d = cheb_diff2( S.x );
S.xu = d(:,:,:,1);
S.xv = d(:,:,:,2);
d = cheb_diff2( S.xu );
S.xuu = d(:,:,:,1);
S.xuv = d(:,:,:,2);
d = cheb_diff2( S.xv );
S.xvv = d(:,:,:,2);

%% initial iterate
tuv0 = zeros(3,1);
if 0
    fid = fopen('tuv.dat','r');
    tuv0(1) = str2num(fgetl(fid));
    tuv0(2:3) = str2num(fgetl(fid))';
    fclose(fid);
end

% perturbation
tuv = [0;0;0];
% tuv = min( 1.0, max( -1.0, tuv0 + 0.1 * ( 2*rand(3,1) - 1) ) );

%% parameters
EPSxyz = 1e-13;
EPSuv = 2.0 * eps('double');

TOLuv = 1e-13;
upperb = ones(3,1) + TOLuv;
lowerb = -upperb;

nitmax = 30;

%% Newton
converged = 0;

res = [];
cnd = [];
for it = 0:nitmax
    % compute residual
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    r = c - s;
    %     r = s - c;
    
    ct = chebval1( C.xt, tuv(1) )';
    su = ICT2unstr( S.xu, tuv(2:3)' )';
    sv = ICT2unstr( S.xv, tuv(2:3)' )';
%     n = cross( su, sv );
    
    %     f = [ dot( n, ct ) ; [su,sv]'*r ]
    f = [ dot( ct, r ) ; -dot( su, r ); -dot( sv, r ) ];
% f = r;
    
    resf = norm( f );
    
    
    % Jacobian matrix
    ctt = chebval1( C.xtt, tuv(1) )';
    suu = ICT2unstr( S.xuu, tuv(2:3)' )';
    suv = ICT2unstr( S.xuv, tuv(2:3)' )';
    svv = ICT2unstr( S.xvv, tuv(2:3)' )';
    
    %     jac = [ ...
    %         dot( n, ctt ) , ...
    %         dot( cross( suu, sv ) + cross( su, suv ), ct ), ...
    %         dot( cross( suv, sv ) + cross( su, svv ), ct ) ;
    %         -[su,sv]' * ct, [suu,suv]' * r, [suv,svv]' * r ]


    hess(1,1) = dot( ctt, r ) + dot( ct, ct );
    hess(1,2) = -dot( ct, su );
    hess(1,3) = -dot( ct, sv );
    hess(2:3,1) = hess(1,2:3);
    hess(2,2) = dot( su, su ) - dot( suu, r );
    hess(2,3) = dot( su, sv ) - dot( suv, r );
    hess(3,3) = dot( sv, sv ) - dot( svv, r );
    hess(3,2) = hess(2,3);


    jac = [ ct, -su, -sv ];
    
%         condJ = cond(jac);
    condJ = cond(hess);
    cnd = [ cnd; condJ ];
    
    % solve for Newton step
%         dtuv = LSsolveQR( jac, -f );
%         dtuv = LSsolveTSVD( jac, -f, TOLuv );
%     dtuv = LSsolveQR( hess, -f );
    dtuv = LSsolveTSVD( hess, -f, TOLuv );
    
    resx = norm( dtuv );
    
    % scale down the step to stay inside the feasible domain
    dtuv = nd_box_reflexions( tuv, lowerb, upperb, dtuv );
    
    % update solution
    tuv = tuv + dtuv;
    
    res = [ res ; [resx, resf, norm(r)]; ];
end

n = size(res,1);

figure;
hold on
plot( 0:n-1, res(:,1), '+-' );
plot( 0:n-1, res(:,2), '.-' );
plot( 0:n-1, res(:,3), 'o-' );
plot( 0:n-1, EPSuv * cnd, '--' );


set( gca, 'yscale', 'log', 'ytick', 10.^(-16:2:4) );%, 'xtick', XTICK );
ylim([1e-16,1]);

legend( {'$|x_{k+1} - x_{k}|$' ; ...
    '$|f(x_k)|$'; ...
    '$|c - s|$'; ...
    '$\epsilon_{x} \, \mathrm{cond}\left(\mathbf{J}_f(x_k)\right)$'}, ...
    'interpreter', 'latex', 'fontsize', 12, ...
    'location', 'northeastoutside' )
grid on

xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 12 )





%%
n = 200;
t = linspace(-1,1,n)';


figure;
hold on

surf_chebyshev( S.x, 1 );

g = chebval1( C.x, t );
plot3( g(:,1), g(:,2), g(:,3), 'k-' );

csol = chebval1( C.x, tuv(1) );
ssol = ICT2unstr( S.x, tuv(2:3)' );

plot3( csol(1), csol(2), csol(3), 'b*' );
plot3( ssol(1), ssol(2), ssol(3), 'r*' );


su = ICT2unstr( S.xu, tuv(2:3)' );
sv = ICT2unstr( S.xv, tuv(2:3)' );
sn = cross( su, sv );
sn = sn / norm( sn );

quiver3( ssol(1), ssol(2), ssol(3), sn(1), sn(2), sn(3), 'r' );
plot3( [csol(1),ssol(1)], [csol(2),ssol(2)], [csol(3),ssol(3)], 'g--' );



axis image vis3d
view(3)
camproj( 'persp' );
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );