clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Newton/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');


tbox = [-1.0000000000000000        1.0000000000000000];
uvbox = [-0.10777282714843750      -0.10775756835937500       -4.9316406250000000E-002  -4.9301147460937500E-002];

%% curve and surface
% C.x = readCoeffs1('../dev_intersection/debugbsi_c.cheb');
C.x = readCoeffs1( 'curv.cheb' );
% C.x = chgvar1(C.x, tbox);
C.xt = cheb_diff1( C.x );
C.xtt = cheb_diff1( C.xt );

% S.x = readCoeffs2('../dev_intersection/debugbsi_s.cheb');
S.x = readCoeffs2( 'surf.cheb' );
% S.x = chgvar2(S.x, reshape(uvbox,2,2));
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
% tuv = [0.5*sum(tbox), 0.5*(uvbox([1,3]) + uvbox([2,4]))]';

% tuv = [0;0;0];
tuv = [0;  -0.107765197753906;  -0.049308776855469];
% tuv = min( 1.0, max( -1.0, tuv0 + 0.1 * ( 2*rand(3,1) - 1) ) );

%% parameters
EPSxyz = 1e-13;
EPSuv = 1e-12;

upperb = ones(3,1) + EPSuv;
lowerb = -upperb;

nitmax = 40;
%% Newton 1
res = [];
for it = 0:nitmax
    % compute residual
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    r = s - c;
    
    % Jacobian matrix
    J(1:3,1) = -chebval1( C.xt, tuv(1) )';
    J(1:3,2) = ICT2unstr( S.xu, tuv(2:3)' )';
    J(1:3,3) = ICT2unstr( S.xv, tuv(2:3)' )';
    ct = -J(1:3,1);
    n = cross(J(1:3,2), J(1:3,3));

    condJ = cond(J);
    
    %     if condJ > EPSuv / eps('double')
    %         % switch to distance minimization with shifted curve
    %         break
    %     end
    
    % solve for Newton step
    dtuv = LSsolveTSVD( J, -r, 3*eps('double') );
    errtuv = norm(dtuv);
    
    % scale down the step to stay inside the feasible domain
    dtuv = nd_box_reflexions( tuv, lowerb, upperb, dtuv );
    
    % update solution
    tuv = tuv + dtuv;
    
    res = [ res ; [errtuv, norm(r), 3*eps('double')*condJ, abs(dot(n,ct))/(norm(ct)*norm(n))]; ];
    
    % termination criteria
    if errtuv < max( EPSuv, 3*eps('double')*condJ )
        % switch to distance minimization with shifted curve
        tuvtmp = tuv;
        break
%         if resx < EPSxyz
%             xyz = 0.5 * ( c + s );
%             converged = 1;
%             break
%         else
%             fprintf('Stagnation...\n');
%         end
    end
end


%% find minimum distance on shifted curve
% ctt = chebval1( C.xtt, tuv(1) )';
% n = cross(J(1:3,2), J(1:3,3));

c = chebval1( C.x, tuv(1) )';
s = ICT2unstr( S.x, tuv(2:3)' )';
norm(c-s)

shiftd = 1e-1;
% sgn = sign(dot(ctt,n));
% shiftv =  sgn * shiftd * n';
% C.x(1,:) = C.x(1,:) + shiftv;

% c0 = C.x(1,:);

hess = zeros(3,3);
for it = 0:nitmax
    ctt = chebval1( C.xtt, tuv(1) )';
    su = ICT2unstr( S.xu, tuv(2:3)' )';
    sv = ICT2unstr( S.xv, tuv(2:3)' )';
    n = cross(su,sv);
    
    sgn = sign(dot(ctt,n));
    shiftv =  sgn * shiftd * n;
%     C.x(1,:) = c0 + shiftv;
    
    % residual
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    r = s - c;
    resx = norm(r);
    r = r - shiftv;
    
    % distance gradient
    ct = chebval1( C.xt, tuv(1) )';
%     su = ICT2unstr( S.xu, tuv(2:3)' )';
%     sv = ICT2unstr( S.xv, tuv(2:3)' )';
%     n = cross(su,sv);
    grad = [-ct, su, sv]'*r;
    
    % Hessian matrix
%     ctt = chebval1( C.xtt, tuv(1) )';
    suu = ICT2unstr( S.xuu, tuv(2:3)' )';
    suv = ICT2unstr( S.xuv, tuv(2:3)' )';
    svv = ICT2unstr( S.xuv, tuv(2:3)' )';
    hess(1,1) = -dot(ctt,r) + dot(ct,ct);
    hess(1,2) = -dot(ct,su);
    hess(1,3) = -dot(ct,sv);
    hess(2,2) = dot(suu,r) + dot(su,su);
    hess(2,3) = dot(suv,r) + dot(su,sv);
    hess(3,3) = dot(svv,r) + dot(sv,sv);
    hess(2:3,1) = hess(1,2:3)';
    hess(3,2) = hess(2,3);
    condJ = cond(hess);
    
    % solve for Newton step
    dtuv = LSsolveTSVD( hess, -grad, 3*eps('double') );
    errtuv = norm(dtuv);
    
    % scale down the step to stay inside the feasible domain
    dtuv = nd_box_reflexions( tuv, lowerb, upperb, dtuv );
    
    % update solution
    tuv = tuv + dtuv;
    
    %     res = [ res ; [errtuv, norm(grad), 3*eps('double')*condJ]; ];
    res = [ res ; [errtuv, resx, 3*eps('double')*condJ, abs(dot(n,ct))/(norm(ct)*norm(n))]; ];
    
    % termination criteria
    if errtuv < max( EPSuv, 3*eps('double')*condJ )
        if resx < EPSxyz
            xyz = 0.5 * ( c + s );
            converged = 1;
            break
        else
            fprintf('Stagnation...\n');
        end
    end
end


n = size(res,1);

if n > 0
    figure;
    hold on
    plot( 0:n-1, res(:,1), '+-' );
    plot( 0:n-1, res(:,2), '.-' );
    plot( 0:n-1, res(:,3), '^--' );
    plot( 0:n-1, res(:,4), 's--' );
    line( [0,n-1], eps('double')*[1,1], 'linestyle', '-' );
    
    set( gca, 'yscale', 'log', 'ytick', 10.^(-16:2:4) );%, 'xtick', XTICK );
    ylim([1e-16,1]);
    
    legend( ...
        {'$|x_{k+1} - x_{k}|$' ; ...
        '$|\mathbf{c} - \mathbf{s}|$'; ...
        '$\epsilon_{x} \, \mathrm{cond}\left(\mathbf{J}_f(x_k)\right)$'; ...
        '$|\mathbf{c}_t \cdot \mathbf{n}|$'}, ...
        'interpreter', 'latex', 'fontsize', 12, ...
        'location', 'northeastoutside' )
    grid on
    
    xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 12 )
end


%%
n = 200;
t = linspace(-1,1,n)';

% C.x(1,:) = c0;

c = chebval1( C.x, tuv(1) )';
s = ICT2unstr( S.x, tuv(2:3)' )';
norm(c-s)



figure;
hold on

% quiver3(c(1), c(2), c(3), T(1), T(2), T(3), 'r');
% quiver3(c(1), c(2), c(3), N(1), N(2), N(3), 'g');
% quiver3(c(1), c(2), c(3), B(1), B(2), B(3), 'b');

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


% plot3(o(1), o(2), o(3), 'kx');

axis image vis3d
view(3)
camproj( 'persp' );
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );