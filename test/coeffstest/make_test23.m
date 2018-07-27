clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 16 );
EPS = 10*eps('double');

rc = 0.4;
angc = 1*pi();
angsu = 0.8*pi();
angsv = 0.4*pi();


%%
cls = linspecerp(2);
nisouv = 5;


%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    r = 0.2*(2.0*rand(2,1) - 1.0);
    switch isurf
        case 1
            [u,v] = tpgrid( 0.5*angsu*cglpoints(-1,1,M), 0.5*angsv*cglpoints(-1,1,M) );
%             u = u + r(1);
            v = v + r(2);
            x = sin(u).*cos(v);
            y = sin(v);
            z = -1 + cos(u).*cos(v);
        case 2
            [u,v] = tpgrid( 0.5*angc*cglpoints(-1,1,M), cglpoints(-1,1,N) );
            u = u + r(1);
            v = v + r(2);
            x = rc*(-1 + cos(u));
            y = rc*sin(u);
            z = v;
    end
    
    c = cat(3, fcht2d(x), fcht2d(y), fcht2d(z));
    c(abs(c) < EPS) = 0;
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,23) );
    
    
    D = cheb_diff2( c );
    x = cat( 3, ifcht2d(c(:,:,1)), ifcht2d(c(:,:,2)), ifcht2d(c(:,:,3)) );
    xu = cat( 3, ifcht2d(D(:,:,1,1)), ifcht2d(D(:,:,2,1)), ifcht2d(D(:,:,3,1)) );
    xv = cat( 3, ifcht2d(D(:,:,1,2)), ifcht2d(D(:,:,2,2)), ifcht2d(D(:,:,3,2)) );
    xn = zeros(size(xu));
    for k = 1:3
        xn(:,:,k) = ...
            xu(:,:,1+mod(k,3)) .* xv(:,:,1+mod(k+1,3)) - ...
            xu(:,:,1+mod(k+1,3)) .* xv(:,:,1+mod(k,3));
    end
    
%     quiver3( x(:,:,1), x(:,:,2), x(:,:,3), xu(:,:,1), xu(:,:,2), xu(:,:,3), 'r' );
%     quiver3( x(:,:,1), x(:,:,2), x(:,:,3), xv(:,:,1), xv(:,:,2), xv(:,:,3), 'g' );
%     quiver3( x(:,:,1), x(:,:,2), x(:,:,3), xn(:,:,1), xn(:,:,2), xn(:,:,3), 'b' );
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


% make_PN_coeff_matrices