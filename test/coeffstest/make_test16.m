clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');

[M,N] = deal( 16 );

R = 0.2;
ang1 = 0.5*pi();
ang2 = 0.3*pi();

shift = R * [0.02,0.0,0];

[u,v] = tpgrid( 0.5 * ang1 * (cglpoints(-1,1,M)+1), 0.5 * cglpoints(-1,1,N) );


%%
u1 = -1;
u2 = 1;
% v12 = 2 * [1.0, -cos(ang2);0, sin(ang2)] \ ...
% [shift(1) - R*sin(ang2)*cos(ang1) ; R*(1 - cos(ang1)*cos(ang2))]
v12 = [-1, cos(ang2); 0,sin(ang2)] \ [ 0.5*cos(ang2)+shift(1) - 0.5 ; 0.5*sin(ang2)];
v12 = 2*v12 - 1;

uv1 = [u1, v12(1)];
uv2 = [u2, v12(2)];

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
            
            uv = uv1;%[-1,0];
        case 2
            x1 = v;
            y1 = R * cos( u+pi()-ang1 );
            z1 = R * sin( u+pi()-ang1 );
%             y1 = R * cos( pi()-u );
%             z1 = R * sin( pi()-u );
            
            x = cos(ang2)*x1 + -sin(ang2)*y1 + shift(1);
            y = sin(ang2)*x1 + cos(ang2)*y1 + shift(2);
            z = z1 + shift(3);
            
            uv = uv2;%[1,0];
    end
    
    C = cat( 3, fcht2d(x), fcht2d(y), fcht2d(z) );
    xi(isurf,:) = ICT2unstr( C, uv );
    
    [si,isou,isov] = surf_chebyshev( C, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( C, sprintf('C%d_test%2.2d.txt',isurf,16) );
    
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
%     quiver3( x(:,:,1), x(:,:,2), x(:,:,3), xn(:,:,1), xn(:,:,2), xn(:,:,3), 'b' );
plot3( xi(isurf,1), xi(isurf,2), xi(isurf,3), '*', 'color', 0.25*cls(isurf,:) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');

  
norm( xi(1,:) - xi(2,:) )

make_PN_coeff_matrices





% pi = [ 0.5, 0, R; 0.5*cos(ang2) + shift(1), 0.5*sin(ang2), R];
% 
% for isurf = 1:2
% plot3( pi(isurf,1), pi(isurf,2), pi(isurf,3), '*', 'color', 0.25*cls(isurf,:) );
% end
