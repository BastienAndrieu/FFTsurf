clc; clear; close all

nisouv = 7;

%%
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/LIBS/mylib/TopologicalEntitites/tests/');
addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/');

%% "Skeleton" surface
%C = readCoeffs2('/stck/bandrieu/Bureau/FFTsurf/cases/Jouke2/reference/Patchs/C_12.cheb');
C = readCoeffs2('c_skeleton.cheb');
B = OBB_Chebyshev2(C);

C(1,1,:), B.extents

fid = fopen('transformation.dat','w');
fprintf(fid, '%d %d %d %d\n', C(1,1,1), C(1,1,2), C(1,1,3), norm(B.extents));
fclose(fid);

C(1,1,:) = 0;
C = C/norm(B.extents);

writeCoeffs2( C, 'c_skeleton.cheb' );


