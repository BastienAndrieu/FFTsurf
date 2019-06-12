clc; clear; close all

x0 = importdata('debug_rip_xyz0.dat');
x1 = importdata('debug_rip_xyz1.dat');
u0 = importdata('debug_rip_uv0.dat');
u1 = importdata('debug_rip_uv1.dat');

%% xyz
figure;
hold on

plot3(x0(:,1), x0(:,2), x0(:,3), 'k-');
plot3(x1(:,1), x1(:,2), x1(:,3), 'r.-');

axis image%equal %vis3d
view(3)

%% uv
figure;
for i = 1:2
    subplot(1,2,i);
    hold on
    plot(u0(:,2*i-1), u0(:,2*i), 'k-');
    plot(u1(:,2*i-1), u1(:,2*i), 'r.-');
    axis equal
end


