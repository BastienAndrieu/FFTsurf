clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/FourierContinuation/Chebyshev/Bezier/');

fid = fopen( 'b.dat', 'r' );
% fprintf( fid, '%d %d\n', size(b,1), size(b,2) );
nd = str2num( fgetl(fid) );
b = zeros(nd(1),nd(2));
for j = 1:nd(2)
    for i = 1:nd(1)
        b(i,j) = str2num( fgetl(fid) );
    end
end
fclose(fid);


bl = importdata('bl.dat');
br = importdata('br.dat');
f = importdata('f.dat');


cl = ['k','r','b'];

figure;
hold on

for i = 1:3
    switch i
        case 1
            bi = b;
        case 2
            bi = bl;
        case 3
            bi = br;
    end
    
    plot( bi(:,1), bi(:,2), '.:', 'color', cl(i) );
    
    g = evalBernstein( bi, linspace(0,1,100)' );
    plot( g(:,1), g(:,2), '-', 'color', cl(i) );
    
    
end

% plot( b(:,1), b(:,2), 'k.-' );
% plot( bl(:,1), bl(:,2), 'r.-' );
% plot( br(:,1), br(:,2), 'b.-' );
plot( f(1), f(2), 'k*' );
axis image
