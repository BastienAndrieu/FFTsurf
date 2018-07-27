clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Newton/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');

%% curve and surface
% C.x = readCoeffs1('../dev_intersection/debugbsi_c.cheb');
C.x = readCoeffs1( 'curv.cheb' );
C.xt = cheb_diff1( C.x );
C.xtt = cheb_diff1( C.xt );

% S.x = readCoeffs2('../dev_intersection/debugbsi_s.cheb');
S.x = readCoeffs2( 'surf.cheb' );
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
EPSuv = 1e-12;

upperb = ones(3,1) + EPSuv;
lowerb = -upperb;

nitmax = 40;

%% Newton
converged = 0;

res = [];
cnd = [];
J = zeros(3,3);

for it = 0:nitmax
    % compute residual
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    r = s - c;
    resx = norm(r);
    
    % Jacobian matrix
    J(1:3,1) = -chebval1( C.xt, tuv(1) )';
    J(1:3,2) = ICT2unstr( S.xu, tuv(2:3)' )';
    J(1:3,3) = ICT2unstr( S.xv, tuv(2:3)' )';
    ct = -J(1:3,1);
    su = J(1:3,2);
    sv = J(1:3,3);
    n = cross(su, sv);
    dotnct = -dot(n,ct);
    
    [U,W,V] = svd(J(1:3,1:3));
    w = diag(W);
    wmax = max(w);
    [wmin,imin] = min(w);
%     wmin / (wmax * eps('double')/EPSuv)
    if wmin < wmax * eps('double')/EPSuv
        ctt = chebval1( C.xtt, tuv(1)   )';
        suu = ICT2unstr(S.xuu, tuv(2:3)')';
        suv = ICT2unstr(S.xuv, tuv(2:3)')';
        svv = ICT2unstr(S.xvv, tuv(2:3)')';

        J(4,:) = [
            dot(ctt, n), ...
            dot(ct, cross(suu, sv) + cross(su, suv)), ...
            dot(ct, cross(suv, sv) + cross(su, svv))
            ];
        r(4) = dotnct;
    else
        J = J(1:3,1:3);
        r = r(1:3);
    end
    
    condJ = cond(J);
    cnd = [ cnd ; condJ ];
    
    abc = [su, sv, n] \ ct;
    
    % solve for Newton step
    dtuv = LSsolveTSVD( J, -r, 3*eps('double') );
    % dtuv = LSsolveTSVD( J, -[r;dotnct], 3*eps('double') );
    errtuv = norm(dtuv);
    
    % scale down the step to stay inside the feasible domain
    dtuv = nd_box_reflexions( tuv, lowerb, upperb, dtuv );
    
    % update solution
    tuv = tuv + dtuv;
    
    res = [ res ; [errtuv, resx, 3*eps('double')*condJ, abs(abc(3))]; ];
    %     res = [ res ; [norm(tuv-tuv0), resx, 3*eps('double')*condJ]; ];
    
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
if 0
    J = zeros(1,3);
    res = [];
    for it = 1:10
        % compute residual
        c = chebval1( C.x, tuv(1) )';
        s = ICT2unstr( S.x, tuv(2:3)' )';
        
        ct = chebval1( C.xt, tuv(1) )';
        su = ICT2unstr( S.xu, tuv(2:3)' )';
        sv = ICT2unstr( S.xv, tuv(2:3)' )';
        n = cross(su, sv);
        
        r = dot(n, ct);
        
        % Jacobian matrix
        ctt = chebval1(C.xtt, tuv(1))';
        suu = ICT2unstr(S.xuu, tuv(2:3)')';
        suv = ICT2unstr(S.xuv, tuv(2:3)')';
        svv = ICT2unstr(S.xvv, tuv(2:3)')';
        
        J = [
            dot(ctt, n), ...
            dot(ct, cross(suu, sv) + cross(su, suv)), ...
            dot(ct, cross(suv, sv) + cross(su, svv))
            ];
        dtuv = -J\r;
        
        res = [ res ; [norm(dtuv), norm(c-s), 3*eps('double')*condJ, abs(r)]; ];
    end
    
    n = size(res,1);
    
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

%% projections
if 0
    c = chebval1(C.x, tuv(1))';
    ct = chebval1(C.xt, tuv(1))';
    ctt = chebval1(C.xtt, tuv(1))';
    o = c + ctt/norm(ctt);
    
    T = ct / norm(ct);
    N = ctt / norm(ctt);
    B = cross(T,N);
    
    su = ICT2unstr(S.xu, tuv(2:3)')';
    sv = ICT2unstr(S.xv, tuv(2:3)')';
    susu = dot(su,su);
    susv = dot(su,sv);
    svsv = dot(sv,sv);
    
    n = cross(su,sv);
    n = n / norm(n);
    q = o - s;
    norm((eye(3) - n*n')*q)
    duv = [susu, susv; susv, svsv] \ ([su,sv]'*(o - s));
    norm(duv)
end

if 0
    nit = 5;
    nits = 1;
    nitc = nits;
    for it = 1:nit
        % project on surface
        c = chebval1( C.x, tuv(1) )';
        for its = 1:nits
            s = ICT2unstr( S.x, tuv(2:3)' )';
            r = c - s;
            resx = norm(r);
            su = ICT2unstr( S.xu, tuv(2:3)' )';
            sv = ICT2unstr( S.xv, tuv(2:3)' )';
            susu = dot(su,su);
            susv = dot(su,sv);
            svsv = dot(sv,sv);
            J = [susu, susv; susv, svsv];
            duv = J \ [dot(su,r); dot(sv,r)];
            erruv = norm(duv);
            %         erruv = norm(tuv-tuv0);
            res = [ res ; [erruv, resx, 2*eps('double')*cond(J)]; ];
            tuv(2:3) = tuv(2:3) + duv;
            n = n + 1;
        end
        
        s = ICT2unstr( S.x, tuv(2:3)' )';
        for itc = 1:nitc
            c = chebval1( C.x, tuv(1) )';
            r = s - c;
            resx = norm(r);
            ct = chebval1( C.xt, tuv(1) )';
            dt = dot(r, ct);
            erruv = abs(dt);
            %                 erruv = norm(tuv-tuv0);
            res = [ res ; [erruv, resx, 2*eps('double')*cond(J)]; ];
            tuv(1) = tuv(1) + dt;
            n = n + 1;
        end
    end
    
    figure;
    hold on
    plot( 0:n-1, res(:,1), '+-' );
    plot( 0:n-1, res(:,2), '.-' );
    % plot( 0:n-1, res(:,3), '--' );
    set( gca, 'yscale', 'log', 'ytick', 10.^(-16:2:4) );%, 'xtick', XTICK );
    ylim([1e-16,1]);
    legend( {'$|x_{k+1} - x_{k}|$' ; ...
        '$|c - s|$'; ...
        '$\epsilon_{x} \, \mathrm{cond}\left(\mathbf{J}_f(x_k)\right)$'}, ...
        'interpreter', 'latex', 'fontsize', 12, ...
        'location', 'northeastoutside' )
    grid on
    xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 12 )
end

%%
n = 200;
t = linspace(-1,1,n)';


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