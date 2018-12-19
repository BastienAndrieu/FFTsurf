clc; clear; close all

u = importdata('bpts.dat');
e = importdata('bedg.dat');
c = readCoeffs2('c.cheb');

x = ICT2unstr(c, u);

l = zeros(size(u,1),2);

figure;
hold on
for i = 1:size(e,1)
    plot(u(e(i,:),1), u(e(i,:),2), 'b.-');
    l(i,1) = norm(u(e(i,1),:) - u(e(i,2),:));
end
axis image


figure;
hold on
surf_chebyshev(c, 1);
for i = 1:size(e,1)
    plot3(x(e(i,:),1), x(e(i,:),2), x(e(i,:),3), 'b.-');
    l(i,2) = norm(x(e(i,1),:) - x(e(i,2),:));
end
axis image vis3d
view(3)
camlight(30,30);


figure;
semilogy(l,'.-');
legend({'u';'x'});