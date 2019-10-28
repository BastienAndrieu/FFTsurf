clc; clear; close all

pth = '../cases/jouke/output/debug/brep/';


verts = importdata([pth, 'verts.dat']);
edges = importdata([pth, 'edges.dat']); % face, orig, prev(2), next(2)

v = [3,8];

xold = importdata('path_old.dat');
xnew = importdata('path_new.dat');
xhe = importdata('hyperedge.dat');

xold = [xold; xold(1,:)];
xnew = [xnew; xnew(1,:)];


% i = 2;
% x = [0.49999999999999956       -4.6259292692714853E-018   0.0000000000000000];
% d = [2.1494334547362480E-015   8.9430006720017195E-003  -0.0000000000000000];

% i = 55;
% x = [6.6819617312574520E-002  0.49551502373005890        0.0000000000000000];
% d = [-1.3587153371560256E-002   1.8322116286615296E-003  -0.0000000000000000];

% y = [0.49992002273898117        8.9426430465015805E-003   0.0000000000000000];

figure;
hold on

plot3(xold(:,1), xold(:,2), xold(:,3), 'k.-');
%plot3(xhe(:,1), xhe(:,2), xhe(:,3), 'b.-');
% plot3(xhe(1,1), xhe(1,2), xhe(1,3), 'go');

% plot3(xold(i,1), xold(i,2), xold(i,3), 'ko');
% quiver3(x(1), x(2), x(3), d(1), d(2), d(3), 0, 'r');
%
% plot3(y(1), y(2), y(3), 'g*');

plot3(xnew(:,1), xnew(:,2), xnew(:,3), 'r.-');

% x = [2.7600000007077942E-002  0.17544229821261245       0.97430446801050541];
% d = [-0.0000000000000000        0.0000000000000000        6.8429410829924806E-002];
%
% quiver3(x(1), x(2), x(3), d(1), d(2), d(3), 0, 'r', 'maxheadsize', 2);

% plot3(verts(v,1), verts(v,2), verts(v,3), 'r*');

m = min(size(xold,1), size(xnew,1)) - 1;
quiver3(xold(1:m,1), xold(1:m,2), xold(1:m,3), ...
    xnew(1:m,1) - xold(1:m,1), xnew(1:m,2) - xold(1:m,2), xnew(1:m,3) - xold(1:m,3), ...
    0, 'g' );

axis equal% vis3d
view(3)
camproj('persp');

xlabel('x');
ylabel('y');
zlabel('z');


