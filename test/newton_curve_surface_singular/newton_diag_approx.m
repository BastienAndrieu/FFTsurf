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
C.xtt = cheb_diff1( C.x );

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
if 1
    fid = fopen('tuv.dat','r');
    tuv0(1) = str2num(fgetl(fid));
    tuv0(2:3) = str2num(fgetl(fid))';
    fclose(fid);
end

% perturbation
tuv = tuv0 + 1e-1*(2*rand(3,1) - 1.0);
% tuv = [0;0;0];
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

stopnow = 0;
for it = 0:nitmax
    if stopnow 
%         xyz = 0.5 * ( c + s );
%             converged = 1;
%             break
    end
    
    % compute residual
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    r = s - c;
    resx = norm(r);
    
    % Jacobian matrix
    J(:,1) = -chebval1( C.xt, tuv(1) )';
    J(:,2) = ICT2unstr( S.xu, tuv(2:3)' )';
    J(:,3) = ICT2unstr( S.xv, tuv(2:3)' )';
    condJ = cond(J);
    cnd = [ cnd ; condJ ];
    
    
    abc = [J(:,2:3), cross(J(:,2), J(:,3))] \ (-J(:,1));
    
    
    % solve for Newton step
    if 0%it > 1 && condJ > EPSuv/eps('double')
        %
    else
        dtuv = LSsolveTSVD( J, -r, 3*eps('double') );
    end
    
    if it > 1 && condJ > EPSuv/eps('double') && ~stopnow
        dtuv = 2*dtuv;
        stopnow = 1;
    end
    
    errtuv = norm(dtuv);
    
    
    % scale down the step to stay inside the feasible domain
    dtuv = nd_box_reflexions( tuv, lowerb, upperb, dtuv );
    
   
    
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
    
    % update solution
    tuv = tuv + dtuv;
    
    rprev = r;
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
    
    legend( {'$|x_{k+1} - x_{k}|$' ; ...
        '$|\mathbf{c} - \mathbf{s}|$'; ...
        '$\epsilon_{x} \, \mathrm{cond}\left(\mathbf{J}_f(x_k)\right)$'; ...
        '$|\mathbf{c}_t \cdot \mathbf{n}|$'}, ...
        'interpreter', 'latex', 'fontsize', 12, ...
        'location', 'northeastoutside' )
    grid on
    
    xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 12 )
end



%%
a = abc(1);
b = abc(2);
res = [];
for it = 1:10
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    
    ct = chebval1( C.xt, tuv(1)   )';
    su = ICT2unstr(S.xu, tuv(2:3)')';
    sv = ICT2unstr(S.xv, tuv(2:3)')';
    
    r = a*su + b*sv - ct;

    
    J(:,1) = -chebval1(C.xtt, tuv(1))';
    suu = ICT2unstr(S.xuu, tuv(2:3)')';
    suv = ICT2unstr(S.xuv, tuv(2:3)')';
    svv = ICT2unstr(S.xvv, tuv(2:3)')';
    J(:,2) = a*suu + b*suv;
    J(:,3) = a*suv + b*svv;
    
    dtuv = -J\r;
    
    res = [ res ; [norm(dtuv), norm(c-s), eps('double')*cond(J), norm(r)]; ];
end


figure;
    hold on
    plot( res(:,1), '+-' );
    plot( res(:,2), '.-' );
    plot( res(:,3), '^--' );
    plot( res(:,4), 's--' );
    line( [0,n-1], eps('double')*[1,1], 'linestyle', '-' );
    
    set( gca, 'yscale', 'log', 'ytick', 10.^(-16:2:4) );%, 'xtick', XTICK );
    ylim([1e-16,1]);
    
    legend( {'$|x_{k+1} - x_{k}|$' ; ...
        '$|\mathbf{c} - \mathbf{s}|$'; ...
        '$|\mathbf{c}_t - a\mathbf{s}_u - b\mathbf{s}_v|$'; ...
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