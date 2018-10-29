clc; clear; close all

addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

nsurf = 6+12+8;
ninstants = 100;
n = 8;

cls = colorcet( 'I2', 'N', nsurf );
cls = cls(randperm(nsurf),:);
cls = CC( cls, 0.0, 0.8, 1.5 );

% axlim = repmat([-1,1],1,3);
% axlim(3:4) = [-0.2,1.2];

R = 0.15;
ctr = [0.5,0.75,0.5];
% axlim = 3.0*R*repmat([-1.0,1.0],1,3);
% for i = 1:3
%     axlim(2*i-1:2*i) = axlim(2*i-1:2*i) + c(i);
% end
axlim = repmat(2.0*[-1,1],1,3);
axlim(1:2) = axlim(1:2) + 0.5;
axlim(3:4) = axlim(3:4) - 0.5;
axlim(5:6) = axlim(5:6) + 0.5;
shift = [-0.05,0,-0.05];
camtg = 0.5*(axlim([1,3,5]) + axlim([2,4,6]));
r_cam = 20.0;
[az,el] = deal(60,30);
[xcam, ycam, zcam] = sph2cart(az*pi()/180.0, el*pi()/180.0, r_cam);
campos = camtg + [xcam, ycam, zcam];


figure('units','normalized','position',[.05 .05 .9 .9]);
for instant = ninstants%:10:ninstants
    %     clf
    hold on
    for i = 1:nsurf
        c = readCoeffs2(sprintf('instant_%3.3dc_%3.3d.cheb',instant,i));
        x = cat(3, ifcht2d(c(:,:,1)), ifcht2d(c(:,:,2)), ifcht2d(c(:,:,3)));
%         for j = 1:3
%             x(:,:,j) = (x(:,:,j) + shift(j) - ctr(j)) / R;
%         end
        surf( x(:,:,1), x(:,:,2), x(:,:,3), ...
            'edgecolor','none','facecolor',cls(i,:),'diffusestrength',0.7, ...
            'facelighting','gouraud','specularstrength',0);
%         si = surf_chebyshev(c, 1, n);
%         set(si, 'facecolor', cls(i,:) , 'specularstrength', 0);
    end
    %     view(30,30);
    %     axis(axlim);
    %     daspect([1,1,1]);
    %     camlight(30,60);
    %     camproj('persp');
    %     grid on
    %     axis vis3d
    %
    %     xlabel('x');
    %     ylabel('y');
    %     zlabel('z');
    %
    %     pause(1e-2);
end
view(30,30);
% axis(axlim);
% daspect([1,1,1]);
axis image
camlight(30,60);
camproj('persp');
grid on
axis vis3d

xlabel('x');
ylabel('y');
zlabel('z');



