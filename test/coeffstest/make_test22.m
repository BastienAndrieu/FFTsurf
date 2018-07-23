clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 16 );
[u,v] = tpgrid( cglpoints(-1,1,M), cglpoints(-1,1,N) );

%%
cls = linspecerp(2);
nisouv = 5;


%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    r = 0.1*(2.0*rand(2,1) - 1.0);
    switch isurf
        case 1
            x = u + r(1);
            y = v + 0;%r(2);
            z = -0.6*(x.^2 + y.^2);
        case 2
%             t = 
            x = 0.5*u + r(1);
            y = 0.7/0.25*x.^2;%0.7*u.^2;
            z = v + r(2);
    end
    
    c = cat(3, fcht2d(x), fcht2d(y), fcht2d(z));
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,22) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


% make_PN_coeff_matrices