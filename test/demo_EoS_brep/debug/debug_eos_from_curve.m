clc; clear; close all

addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

g = importdata('g.dat');

figure; hold on

plot3(g(:,1), g(:,2), g(:,3), 'k.-');

for i = 1:2
    c = readCoeffs2(sprintf('enve_rl_%d.cheb',i));
    surf_chebyshev(c, 1);
    
%     u = importdata(sprintf('uvc_%d.dat',i));
%     x = ICT2unstr(c, u);
%     plot3(x(:,1), x(:,2), x(:,3), 'ro-');
end

y = importdata('fitedge_y.dat');
x = y(:,5:7);
% plot3(x(:,1), x(:,2), x(:,3), 'ro-');

% c = importdata('fitedge_c.dat');
% x = ifcht(c(:,5:7));
% plot3(x(:,1), x(:,2), x(:,3), 'ro');
% t = linspace(-1,1,1000)';
% x = chebval1(c(:,5:7), t);
% plot3(x(:,1), x(:,2), x(:,3), 'r-');

% g = importdata('gc.dat');
% dg = importdata('dgc.dat');
% quiver3(g(:,1), g(:,2), g(:,3), dg(:,1), dg(:,2), dg(:,3), 'r');

for i = 1:2
    e = importdata(sprintf('erlc_%d.dat',i));
    plot3(e(:,1), e(:,2), e(:,3), 'ro-');
end


% for i = 1:3
%     X{i} = importdata(sprintf('xyz_%d.dat',i));
% end
% surf(X{1}, X{2}, X{3}, 'facecolor', 'y', 'specularstrength', 0);
X = readCoeffs2('xyzc.dat');
surf(X(:,:,1), X(:,:,2), X(:,:,3), 'facecolor', 'y', 'specularstrength', 0);

axis image vis3d
% axis(repmat(0.8*[-1,1],1,3));
% daspect([1,1,1])
view(130,15);
camlight(-30,30);
camproj('persp');


% return
% %%
% x = importdata('fitedge_x.dat');
% figure;
% plot(linspace(0,1,length(x)), x,'.-');

%%
figure;
for i = 1:2
    subplot(1,2,i);
    hold on
    
    u = y(:,2*i-1:2*i);
    plot(u(:,1), u(:,2), '.-');
    
    u = importdata(sprintf('uvc_%d.dat',i));
    plot(u(:,1), u(:,2), 'ro-');
    
    axis(repmat(1.05*[-1,1], 1, 2));
    daspect([1,1,1]);
end