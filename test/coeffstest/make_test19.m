clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');


% f = 5;
% x = [
%     -1,-1,-1,-1;
%     1 , 1, 1, 1;
%     -1,-1,-1,-1;
%     1 , 1, 1, 1];
% x(2:3,:) = f*x(2:3,:);
% 
% y = x';
% % y = repmat( linspace(0,1,4), 4, 1 );
% 
% v = linspace(-1,1,4);
% z = 0.05 * ( repmat(v',1,4) + repmat(v,4,1) );
% 
% b = cat( 3, x, y, z );
% 
% 
% n = 200;
% t = linspace( 0, 1, n )';
% 
% figure;
% hold on
% 
% % surf( b(:,:,1), b(:,:,2), b(:,:,3), ...
% %  'edgecolor','r','facecolor','none');
% 
% s = eval_Bezier_patch_tpgrid( b, t, t );
% surf( s(:,:,1), s(:,:,2), s(:,:,3), ...
%     'edgecolor','none','facecolor','y','diffusestrength',0.7, ...
%     'facelighting','gouraud','specularstrength',0.15);
% % contourSurf( s(:,:,1), s(:,:,2), s(:,:,3) );
% 
% 
% 
% 
% xc = [ 0, 0, 0, 0;
%     1, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0];
% yc = [ 0, 1, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0];
% scale2 = 0.9;
% 
% zc = [ 0.0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 0;
%     0, 0, 0, 1];
% xc = scale2*xc;
% yc = scale2*yc;
% 
% 
% 
% c = cat( 3, xc, yc, zc );
% b2 = chebyshev2bezier_2(c)
% 
% [si,isou,isov] = surf_chebyshev( c, 1, n );
% 
% 
% 
% x = repmat( scale2*linspace(-1,1,4)', 1, 4 );
% y = x';
% fz = 0.6;
% z = 0.6 * [
%     1    -5     5    -1
%     -5    25   -25     5
%     5   -25    25    -5
%     -1     5    -5     1];
% 
% 
% 
% 
% 
% axis image vis3d
% view(3)
% camlight(30,30);
% 
% return



%%

scale2 = 0.9;

%%
cls = linspecerp(2);
nisouv = 10;

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            f = 5;
            x = [
                -1,-1,-1,-1;
                1 , 1, 1, 1;
                -1,-1,-1,-1;
                1 , 1, 1, 1];
            x(2:3,:) = f*x(2:3,:);            
            y = x';
            
            v = linspace(-1,1,4);
            z = 0.05 * ( repmat(v',1,4) + repmat(v,4,1) );
            
        case 2
            x = repmat( scale2*linspace(-1,1,4)', 1, 4 );
            y = x';
            fz = 0.6;
            z = 0.6 * [
                1    -5     5    -1
                -5    25   -25     5
                5   -25    25    -5
                -1     5    -5     1];

    end
    
    b = cat( 3, x, y, z );
    c = 0*b;
    A = C2Bmatrix(4);
    At = A';
    for i = 1:3
        c(:,:,i) = A \ b(:,:,i) / At;
    end
    
%     C = cat( 3, xc, yc, zc );
    
    [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
    set( si, 'facecolor', cls(isurf,:) );
    if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,19) );
    
    %     B = chebyshev2bezier_2( C );
    %     surf( B(:,:,1), B(:,:,2), B(:,:,3), 'facecolor', 'none', 'edgecolor', 0.2*cls(isurf,:) );
    
end
axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


make_PN_coeff_matrices