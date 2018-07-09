clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 16 );

R1 = 0.5;
R2 = 1.0;
ang = pi()/6;
H = 0.5;

[u,v] = tpgrid( cglpoints(0,ang,M), cglpoints(0,1,N) );

%%
cls = linspecerp(2);
nisouv = 10;

n = 20;
t = linspace(-1,0,n)';
uv = cat( 3, [t,ones(n,1)], [t,-ones(n,1)] );
mrk = ['.','o'];

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            r = (1-v)*R2 + v*R1;
            x = r .* cos( u );
            y = r .* sin( u );
            z = H*ones(M,N);
        case 2
            x = R2 * cos( u );
            y = R2 * sin( u );
            z = H * v;
    end
    
    C = cat( 3, fcht2d(x), fcht2d(y), fcht2d(z) );
    
    [si,isou,isov] = surf_chebyshev( C, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( C, sprintf('C%d_test%2.2d.txt',isurf,15) );
    
    c = ICT2unstr( C, uv(:,:,isurf) );
    plot3( c(:,1), c(:,2), c(:,3), 'y-', 'marker', mrk(isurf) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
    
make_PN_coeff_matrices