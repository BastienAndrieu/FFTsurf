clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 16 );

R = 0.2;
ang1 = 0.5*pi();
ang2 = 0.3*pi();

shift = R * [0.02,0.0,0];

[u,v] = tpgrid( 0.5 * ang1 * cglpoints(0,2,M), 0.5 * cglpoints(-1,1,N) );

%%
cls = linspecerp(2);
nisouv = 10;


%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            x = v;
            y = R * cos( u );
            z = R * sin( u );
            
            uv = [-1,1];
        case 2
            x1 = v;
            y1 = R * cos( pi()-u );
            z1 = R * sin( pi()-u );
            
            x = cos(ang2)*x1 + -sin(ang2)*y1 + shift(1);
            y = sin(ang2)*x1 + cos(ang2)*y1 + shift(2);
            z = z1 + shift(3);
            
            uv = [1,1];
    end
    
    C = cat( 3, fcht2d(x), fcht2d(y), fcht2d(z) );
    xi(isurf,:) = ICT2unstr( C, uv );
    
    [si,isou,isov] = surf_chebyshev( C, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( C, sprintf('C%d_test%2.2d.txt',isurf,16) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');

  
norm( xi(1,:) - xi(2,:) )

make_PN_coeff_matrices