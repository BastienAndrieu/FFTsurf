clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

shape = 'freeform';% 'sphere';% 'cone';% 'torus';% 

N = 16;

switch shape
    case 'sphere'
        T = 0.7*pi();
        P = 0.3*pi();
                
        u = T*cglpoints(1,-1,N)';
        v = P*cglpoints(1,-1,N);
        
        x = cos(u)*cos(v);
        y = sin(u)*cos(v);
        z = repmat(sin(v), N, 1);
        
    case 'torus'
        R1 = 1.0;
        R2 = 0.4;
        T = 0.9*pi();
        P = 0.9*pi();
        
        u = T*cglpoints(1,-1,N)';
        v = P*cglpoints(1,-1,N);
        
        x = cos(u)*(R1 + R2*cos(v));
        y = sin(u)*(R1 + R2*cos(v));
        z = repmat(R2*sin(v), N, 1);
        
    case 'cone'
        R1 = 1.0;
        R2 = 0.1;
        H = 2.0;
        T = 0.9*pi();
        
        u = T*cglpoints(1,-1,N)';
        v = cglpoints(1,0,N);
        
        r = (1 - v)*R1 + v*R2;
        x = cos(u)*r;
        y = sin(u)*r;
        z = repmat(H*v, N, 1);
        
    case 'freeform'
        %c = readCoeffs2('../demo_EoS_MAT/c_skeleton.cheb');
        c = readCoeffs2('../coeffstest/C1_test09.txt');
        x = ifcht2d(c(:,:,1));
        y = ifcht2d(c(:,:,2));
        z = ifcht2d(c(:,:,3));
end


c = cat(3, fcht2d(x), fcht2d(y), fcht2d(z));

writeCoeffs2(c, 'c.cheb');




figure;
hold on

[S, lu, lv] = surf_chebyshev(c, 1, 200, 9);

axis image vis3d
view(23,16)
camlight(30,30);
camproj('persp');