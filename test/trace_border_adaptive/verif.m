clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

lbshow = 1;
lbscale = 5e-2;
lboffset = 1e-4;
lbuv = -0.9*[1,1];



c = readCoeffs2('surf_x.cheb');

for i = 1:4
    bi = importdata(sprintf('border_%d.dat',i));
    border(i).xyz = bi(:,1:3);
    border(i).uv = bi(:,4:5);
end

% cl = ['r', 'g', 'b', 'c'];
cl = linspecer(4);

figure; 
hold on

surf_chebyshev(c, 1);
if lbshow
    plot_local_basis(c, lbuv, lbscale, lboffset);
end

for i = 1:4
    plot3(border(i).xyz(:,1), border(i).xyz(:,2), border(i).xyz(:,3), '.-', 'color', cl(i,:));
    
%     x = ICT2unstr(c, border(i).uv);
%     plot3(x(:,1), x(:,2), x(:,3), 'o', 'color', cl(i,:));
end

axis image vis3d
view(3)
camproj('persp')
camlight(30,30);


figure; 
hold on

for i = 1:4
    plot(border(i).uv(:,1), border(i).uv(:,2), '.-', 'color', cl(i,:));
end

axis(repmat(1.05*[-1,1],1,2));
daspect([1,1,1])
