clc; clear; close all

c = readCoeffs1('chebfit1_c.cheb');
y = importdata('chebfit1_y.cheb');

n = 1000;
t = linspace(-1,1,n)';

z = chebval1(c, t);



figure;
hold on

plot3(y(:,1), y(:,2), y(:,3), 'b.-');
plot3(z(:,1), z(:,2), z(:,3), 'r-');

view(3);
daspect([1,1,1])
axis vis3d