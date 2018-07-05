clc; clear; %close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Newton/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');

%% curve and surface
C.x = readCoeffs1( 'curv.cheb' );
C.xt = cheb_diff1( C.x );

S.x = readCoeffs2( 'surf.cheb' );
d = cheb_diff2( S.x );
S.xu = d(:,:,:,1);
S.xv = d(:,:,:,2);

%% initial iterate
tuv0 = zeros(3,1);
fid = fopen('tuv.dat','r');
tuv0(1) = str2num(fgetl(fid));
tuv0(2:3) = str2num(fgetl(fid))';
fclose(fid);

% perturbation
tuv = [0;0;0];
% tuv = min( 1.0, max( -1.0, tuv0 + 0.5 * ( 2*rand(3,1) - 1) ) );

%% parameters
EPSxyz = 1e-13;
EPSuv = 2.0 * eps('double');
TOLuv = 1e-10;
seuil_cond = TOLuv / eps('double');

TOLuvbounds = 1e-13;
upperb = ones(3,1) + TOLuvbounds;
lowerb = -upperb;

nitmax = 30;
skipit = 4;
acceleration = 0;

%% Newton
converged = 0;
itconv = nitmax;

muit = zeros(3,1);
qit = zeros(3,1);

condJ = 1;
normdtuvprev = 0;

res = [];
normd = [];
cnd = [];
err = [];
cndf = [];
err_inv = [];


tuvit = [];

for it = 0:nitmax
    tuvit = [ tuvit, tuv ];
    err = [ err; norm( tuv - tuv0 ) ];
    
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    
    r = s - c;
    res = [ res ; norm(r) ];
    
    % Jacobian matrix
    ct = chebval1( C.xt, tuv(1) )';
    su = ICT2unstr( S.xu, tuv(2:3)' )';
    sv = ICT2unstr( S.xv, tuv(2:3)' )';
    
    J = [ -ct, su, sv ];
    condJ = cond(J);
    cnd = [ cnd ; condJ ];
    if 1%condJ < seuil_cond
        Jsafe = J;
        condJsafe = condJ;
    end
    
    % convergence criterion
    if it > 1
        if norm( r ) < EPSxyz && errtuv < condJsafe * EPSuv
            xyz = 0.5 * ( c + s );
            converged = 1;
            itconv = min(it,itconv);
            break
        end
    end
    
    % solve for Newton step
    dtuv = LSsolveQR( Jsafe, -r ); % s_k
    if mod(it,3) == 0
        dtuv = 2*dtuv;
    end
    
    err_inv = [ err_inv; norm( J*dtuv + r ) / ( norm(J) * norm(dtuv) ) ];
    
    normdtuv = norm(dtuv);
    
    errtuv = normdtuv;
    
    normd = [ normd ; normdtuv ];
    
    dtuv = nd_box_reflexions( tuv, lowerb, upperb, dtuv );
    
    % update solution
    tuv = tuv + dtuv; % x_k+1 = x_k + s_k
    
    normdtuvprev = normdtuv;
end

if 0%converged
    for i = 1:size(tuvit,2)
        err(i) = norm( tuvit(:,i) - tuv );
    end
end

n = it+1;

figure('units','normalized','position',[.2,.2,.6,.6] );
hold on
plot( 0:n-1, err, 'o-' );
plot( 0:n-1, res, '.-' );
plot( 0:length(normd)-1, normd, '+-' );
plot( 0:n-1, EPSuv * cnd, '--' );
line( [0,n-1], EPSxyz*[1,1], 'linestyle', '--' );
if 0%converged
    line( itconv*[1,1], [1e-15, 1e1], 'linestyle', '--', 'color', 'k' );
end

xstride = max( 1, round( n / 10 ) );
XTICK = 0:xstride:nitmax;

set( gca, 'yscale', 'log', 'ytick', 10.^(-16:2:4), 'xtick', XTICK );
ylim([1e-16,1]);

legend( {'$|x_k - x_*|$' ; ...
    '$|f(x_k)|$' ; ...
    '$|x_{k+1} - x_{k}|$' ; ...
    '$\epsilon_{x} \, \mathrm{cond}\left(\mathbf{J}_f(x_k)\right)$'; ...
    '$\epsilon_{f}$'}, ...
    'interpreter', 'latex', 'fontsize', 12, ...
    'location', 'northeastoutside' )
grid on

xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 12 )


return
% figure('name','erreur inverse')
% semilogy( err_inv, '.-' );
% grid on

n = length(normd);
q = log( normd(2:n) ) ./ log( normd(1:n-1) );
mu = normd(2:n) ./ normd(1:n-1);

figure('name','order & rate of convergence')
hold on
% plot( 1:n-1, q, 'o-' );
plot( 1:n-1, mu, 'x-' );
if p == 0
    qth = 2;
    muth = 0;
else
    qth = 1;
    muth = p / (p + 1);
end
% plot( [1,n-1], qth*[1,1], 'linestyle', '--' );
plot( [1,n-1], muth*[1,1], 'linestyle', '-.' );
grid on
% legend( ...
%     {'$\frac{\log|s_k|}{\log|s_{k-1}|}$';
%     '$\frac{|s_k|}{|s_{k-1}|}$';
%     '$q$';
%     '$\mu$'}, ...
%     'interpreter', 'latex', 'fontsize', 16, ...
%     'location', 'northeastoutside' )
legend( ...
    {'$\frac{|s_k|}{|s_{k-1}|}$';
    '$\mu$'}, ...
    'interpreter', 'latex', 'fontsize', 16, ...
    'location', 'northeastoutside' )
xlabel( '$k$', 'interpreter', 'latex', 'fontsize', 12 )

if p == 0
    set( gca, 'yscale', 'log' );
end









%% plot 3d
cls = colorcet( 'I2', 'N', 2 );
n = 300;
t = linspace(-1,1,n)';

figure;
hold on

si = surf_chebyshev( S.x, 1, n, 0, 0 );
set( si, 'specularstrength', 0, 'facecolor', cls(1,:) );

g = chebval1( C.x, t );
plot3( g(:,1), g(:,2), g(:,3), 'k-' );

% s0 = ICT2unstr( S.x, tuv0(2:3)' );
% plot3( s0(1), s0(2), s0(3), '.', 'color', 0.5*cls(1,:) );
% g0 = chebval1( C.x, tuv0(1) );
% plot3( g0(1), g0(2), g0(3), 'k.' );

if converged
    plot3( xyz(1), xyz(2), xyz(3), 'k*' );
    axis( xyz([1,1,2,2,3,3])' + repmat( 0.15 * [-1,1], 1, 3 ) );
    daspect([1,1,1]);
else
    axis image
end

axis vis3d
view(3)
camproj( 'persp' );
% camlight(30,30);
[azl,ell] = deal( 120, 60 );
[xl,yl,zl] = sph2cart( pi()*azl/180, pi()*ell/180, 10 );
light( 'style', 'infinite', 'position', [xl,yl,zl] );
light( 'style', 'infinite', 'position', [-xl,-yl,-0.5*zl], 'color', 0.7*[1,1,1] );
