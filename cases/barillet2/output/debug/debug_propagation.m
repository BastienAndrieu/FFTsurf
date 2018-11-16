clc; clear; close all

addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/Matlab/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Chebyshev/');
addpath('/stck8/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

nsurf = 7;
ninstants = 40;
n = 16;

cls = colorcet( 'I2', 'N', nsurf );
cls = cls(randperm(nsurf),:);
cls = CC( cls, 0.0, 0.8, 1.5 );

fa = 0.5*ones(nsurf,1);
% fa([4,5]) = 1;
 fa([1,2]) = 1;
 
axlim = [0.6*repmat([-1,1],1,2), 0, 1];
% axlim = [0,1.25,-0.5,0.5,0,0.75];

figure('units','normalized','position',[.05 .05 .9 .9]);
for instant = 1:1:ninstants
        clf
    hold on
    for i = 1:nsurf
        c = readCoeffs2(sprintf('instant_%3.3dc_%3.3d.cheb',instant,i));
        x = cat(3, ifcht2d(c(:,:,1)), ifcht2d(c(:,:,2)), ifcht2d(c(:,:,3)));
%         for j = 1:3
%             x(:,:,j) = (x(:,:,j) + shift(j) - ctr(j)) / R;
%         end
        surf( x(:,:,1), x(:,:,2), x(:,:,3), ...
            'edgecolor','none','facecolor',cls(i,:),'diffusestrength',0.7, ...
            'facelighting','gouraud','specularstrength',0,'facealpha',fa(i));
%         si = surf_chebyshev(c, 1, n);
%         set(si, 'facecolor', cls(i,:) , 'specularstrength', 0, 'facealpha', fa(i));
    end
        view(30,30);
        axis(axlim);
        daspect([1,1,1]);
        camlight(30,60);
        camproj('persp');
        grid on
        axis vis3d
    
        xlabel('x');
        ylabel('y');
        zlabel('z');
    
        pause(1e-2);
end
% view(30,30);
% % axis(axlim);
% % daspect([1,1,1]);
% axis image
% camlight(30,60);
% camproj('persp');
% grid on
% axis vis3d
% 
% xlabel('x');
% ylabel('y');
% zlabel('z');



