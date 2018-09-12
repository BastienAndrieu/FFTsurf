clc; clear; close all

%%
u = importdata('bpts.dat');
e = importdata('bedg.dat');

figure;
hold on

plot(u(:,1), u(:,2), '.');
for i = 1:size(u,1)
    text(u(i,1), u(i,2), num2str(i), 'color', 'b');
end

for i = 1:size(e,1)
    plot(u(e(i,:),1), u(e(i,:),2), 'r-');
end
axis image



%%
t = importdata('tri.dat');
x = importdata('uv.dat');

figure;
hold on

triplot(t, x(:,1), x(:,2), 'color', 'k');

axis image
