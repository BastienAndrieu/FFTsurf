clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');


x = importdata('eos_from_polyline_x.dat');
y = importdata('eos_from_polyline_y.dat');
z = importdata('eos_from_polyline_z.dat');

err{1} = importdata('eos_from_polyline_err_rad.dat');
err{2} = importdata('eos_from_polyline_err_dot.dat');

for i = 1:2
    figure;
    hold on
    
    surf(x, y, z, max(-17,log10(err{i})), 'facecolor', 'interp', 'edgecolor', 'none');
    
    axis image vis3d
    view(3)
    colormap(parula(17));
    colorbar
    caxis([-16,0]);
end
