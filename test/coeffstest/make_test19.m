clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/Bezier/');
addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

%%
cls = linspecerp(2);
nisouv = 33;
clcnt = ['r','b'];

scale2 = 0.9;


%%
xyz = [
    1.2824729401064705       -1.2824729401064714       -1.2824729401064428E-002
    0.67474999436144178      -0.90000000000000002       0.34812459746238966
    ];
uv = [
    -0.51298917604257710       0.51298917604257699
    0.74972221595715760       -1.0000000000000000
    ];
n = [
    -0.98224937479061902        9.8414524536714593E-002 -0.15968953341922740     
   0.57787383582322160       0.78261675177713919       0.23145809493028410
   ];

%%
figure('units','normalized','position',[.1 .1 .8 .8]);
hold on

for isurf = 1:2
    switch isurf
        case 1
            if 0
            f = 5;
            if 0
                x = [
                    -1,-1,-1,-1;
                    1 , 1, 1, 1;
                    -1,-1,-1,-1;
                    1 , 1, 1, 1];
                x(2:3,:) = f*x(2:3,:);
            else
                x = [
                    -1,-1,-1,-1;
                    5 , 5, 5, 5;
                    -5,-5,-5,-5;
                    1 , 1, 1, 1];
            end
            y = x';
            
            v = linspace(-1,1,4);
            z = 0.08 * ( repmat(v',1,4) + 0.4*repmat(v,4,1) );
            
            else
                
                x = [
                    -1.5    0.4   -3.0   -0.3
                    -1.5    0.8   -2.7   -0.2
                    -1.5    1.9   -4.0   -0.6
                    -1.5    1.0   -2.3   -0.3
                    ];
                
                y = [
                    -1.5   -1.5   -1.5   -1.5
                    1.3    1.3    1.3    1.3
                    -3.5   -3.5   -3.5   -3.5
                    -0.7   -0.7   -0.7   -0.7
                    ];
                
                z = [
                    -0.1   0.0   -0.5   -0.4
                    -0.1   -0.1   -0.4   -0.4
                    0.0   0.0   -0.4   -0.3
                    0.0   0.0   -0.5   -0.3
                    ];
                
                x = x - mean(x(:));
                y = y - mean(y(:));
                z = z - mean(z(:));
                
                x = 2.5*x;
                y = 2.5*y;
                z = 0.5*z;
            end
            
            
        case 2
            x = repmat( scale2*linspace(-1,1,4)', 1, 4 );
            y = x';
            fz = 0.6;
            z = fz * [
                1    -5     5    -1
                -5    25   -25     5
                5   -25    25    -5
                -1     5    -5     1];
            z = z - 0.01;
    end
    
    b = cat( 3, x, y, z );
    if isurf == 1
        b = readCoeffs2('B19.txt');
        bm = mean(mean(b,2),1);
        scl = [2.5,2.5,0.35];
        for idim = 1:3
            b(:,:,idim) = scl(idim) * ( b(:,:,idim) - bm(idim) );
        end
    end
    c = 0*b;
    A = C2Bmatrix(4);
    At = A';
    for i = 1:3
        c(:,:,i) = A \ b(:,:,i) / At;
    end
    
    %     C = cat( 3, xc, yc, zc );
    
    if 0
        [si,isou,isov] = surf_chebyshev( c, 1, 300, nisouv );
        set( si, 'facecolor', cls(isurf,:) );
        if nisouv > 0; set( [isou;isov], 'color', 0.5*cls(isurf,:) ); end
    else
        [si,~,~,~,cnt] = surf_chebyshev( c, 1, 300, 0, 0, 1 );
        set( si, 'facecolor', 'w' );
        set( cnt, 'color', clcnt(isurf) );
    end
    writeCoeffs2( c, sprintf('C%d_test%2.2d.txt',isurf,19) );
    
    %     B = chebyshev2bezier_2( C );
    %     surf( B(:,:,1), B(:,:,2), B(:,:,3), 'facecolor', 'none', 'edgecolor', 0.2*cls(isurf,:) );
    
    if 0% isurf == 1;
        surf( b(:,:,1), b(:,:,2), b(:,:,3), 'facecolor', 'none', 'edgecolor', 0.2*cls(isurf,:) );
    end
    
%     x = ICT2unstr(c, uv(isurf,:));
%     norm(x - xyz(isurf,:))
%     quiver3(x(1), x(2), x(3), n(isurf,1), n(isurf,2), n(isurf,3), 'color', 0.5*cls(isurf,:) );
    
%     pn = readCoeffs2(sprintf('PN%d_test%2.2d.txt',isurf,19) );
%     m = 20;
%     t = linspace(-1,1,m)';
%     x = ICT2tpgrid(c, t, t);
%     q = ICT2tpgrid(pn, t, t);
%     quiver3(x(:,:,1), x(:,:,2), x(:,:,3), q(:,:,1), q(:,:,2), q(:,:,3), 'color', 0.5*cls(isurf,:) );
end


axis image vis3d
view(3)
camproj('persp');
camlight(30,30);
xlabel('x'); ylabel('y'); zlabel('z');


% make_PN_coeff_matrices