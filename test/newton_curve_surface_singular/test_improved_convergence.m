clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Newton/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Experiments/parameters_qT/bivariate/');

C.x = readCoeffs1( 'curv.cheb' );
C.xt = cheb_diff1( C.x );

S.x = readCoeffs2( 'surf.cheb' );
d = cheb_diff2( S.x );
S.xu = d(:,:,:,1);
S.xv = d(:,:,:,2);

tuv0 = zeros(3,1);
fid = fopen('tuv.dat','r');
tuv0(1) = str2num(fgetl(fid));
tuv0(2:3) = str2num(fgetl(fid))';
fclose(fid);
tuv0 = tuv0 + 0.05 * ( 2*rand(3,1) - 1);

%%
tolsv = 1e3 * eps('double');
nitmax = 30;

EPSuv = 1e-13;
EPSxyz = 1e-13;
upperb = ones(3,1) + EPSuv;
lowerb = -upperb;

tuv = tuv0;
res = [];
normd = [];
cnd = [];
for it = 1:nitmax
    c = chebval1( C.x, tuv(1) )';
    s = ICT2unstr( S.x, tuv(2:3)' )';
    
    r = s - c;
    
    res = [ res ; norm(r) ];
    
    
    ct = chebval1( C.xt, tuv(1) )';
    su = ICT2unstr( S.xu, tuv(2:3)' )';
    sv = ICT2unstr( S.xv, tuv(2:3)' )';
    
    J = [ -ct, su, sv ];
    condJ = cond(J);
    cnd = [ cnd ; condJ ];
    
    if norm(r) < EPSxyz
        converged = 1;
        xyz = 0.5 * ( c + s );
        break
    end
    
    %         dtuv = - J \ r;
    if 0
        dtuv = LSsolveQR( J, -r );
    else
        [U,Sv,V] = svd(J);
        Sinv = zeros(size(Sv));
        for i = 1:size(Sv,1)
            if Sv(i,i) > tolsv * Sv(1,1)
                Sinv(i,i) = 1 / Sv(i,i);
            end
        end
        dtuv = - V * Sinv * U' * r;
    end

    e = J * dtuv + r;
    if norm(e) > eps('double') %condJ > EPSxyz/eps('double')
        norm(e)
        correc = - LSsolveQR( J, e );
        dtuv = dtuv + correc;
    end
    
    % scale down Newton step to keep the solution inside feasible region
    lambda = nd_box_constraint( tuv, lowerb, upperb, dtuv );
    dtuv = lambda * dtuv;
    if mod(it,3) == 0
        dtuv = 2.0 * dtuv;
    end
    normd = [ normd ; norm(dtuv) ];
    
    
    
    tuv = tuv + dtuv;
end
n = length(res);

figure
hold on
plot( 1:n, res, '.-' );
plot( 0.5+(1:n-1), normd, '+-' );
plot( 1:n, eps('double')*cnd, '--' );
set( gca, 'yscale', 'log' );
grid on

