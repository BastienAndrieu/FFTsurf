clc; clear; close all


tri = importdata('../connect.dat');
xyz = importdata('../../../../tmp/xyz.dat');
ener = importdata('../../../../tmp/ener.dat');
% grad = importdata('../../../../tmp/grad.dat');
% hess = importdata('../../../../tmp/hess.dat');

l = 3984;%6006;
instant = 1;
nsurf = 6;
cls = colorcet( 'I2', 'N', nsurf );
cls = cls(randperm(nsurf),:);
cls = CC( cls, 0.0, 0.8, 1.5 );


figure('units','normalized', 'position', [.1 .1 .8 .8]);
hold on
% 
for i = 1:6
    if instant == 0
        c = readCoeffs2(sprintf('../../init/coef/c_%3.3d.cheb',i));
    else
        c = readCoeffs2(sprintf('instant_%3.3dc_%3.3d.cheb',instant,i));
    end
    si = surf_chebyshev(c, 1);
    set(si, 'facecolor', cls(i,:) , 'specularstrength', 0);
end


% trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), ener );
trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'none', 'specularstrength', 0);

p = [
    0.48099421103197926       0.85402774834096751       0.39125848138465685     
  0.48599860440472831       0.85455692140582107       0.39085313643468900 
];

d = [5.6532289429916127E-006  -2.6122923885606248E-005  -2.6802247318088726E-005];
q = xyz(l,:) + d;

x2 = [0.48599860440472831       0.85455692140582107       0.39085313643468900];
d2 = [-6.8394673571856124E-005  -1.9682259498634046E-005  -8.4175279719656373E-006];
q2 = x2 + d2;

plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'go');
quiver3(xyz(l,1), xyz(l,2), xyz(l,3), d(1), d(2), d(3), 0, 'r', 'maxheadsize', 10);

quiver3(x2(1), x2(2), x2(3), d2(1), d2(2), d2(3), 0, 'b', 'maxheadsize', 10);

plot3(q(1), q(2), q(3), 'r*');
plot3(q2(1), q2(2), q2(3), 'b*');

plot3(p(:,1), p(:,2), p(:,3), 'yo-');

% axis image
axlim = xyz(l,[1,1,2,2,3,3]);
axlim = axlim + 1e-2*repmat([-1,1],1,3);
axis(axlim);
daspect([1,1,1]);
axis vis3d

view(240,-50)
camproj('persp');
colormap(cbrewer('seq','OrRd',16));
camlight(30,30);

