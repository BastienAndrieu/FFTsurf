clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 4 );
[u,v] = tpgrid( cglpoints(-1,1,M) );

%%
cls = linspecerp(2);
nisouv = 10;

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            x = u + 0.1*v.^2;
            y = v;
            z = -u.^2;% + y.^2;
        case 2
            x = u + 0.1*v.^2;
            y = -v;
            z = u.^2;%-x.^2;
    end
    
    C = cat( 3, fcht2d(x), fcht2d(y), fcht2d(z) );
    
    [si,isou,isov] = surf_chebyshev( C, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( C, sprintf('C%d_test%2.2d.txt',isurf,24) );
    
    D = cheb_diff2( C );
    x = cat( 3, ifcht2d(C(:,:,1)), ifcht2d(C(:,:,2)), ifcht2d(C(:,:,3)) );
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
    quiver3( x(:,:,1), x(:,:,2), x(:,:,3), xn(:,:,1), xn(:,:,2), xn(:,:,3), 'b' );

end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');

  
% make_PN_coeff_matrices

