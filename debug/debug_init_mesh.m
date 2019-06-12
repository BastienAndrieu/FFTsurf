clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/Regression3D/GPeyre/toolbox_graph/');

e = [
    3789           1
    ];
v = [476];

for i = 0:1
    
    tri = importdata(sprintf('connect_%2.2d.dat',i));
    xyz = importdata(sprintf('pos_%3.3d.dat',i));
    
    
    fprintf('%d verts, %d triangles\n',size(xyz,1), size(tri,1));
    
    figure('units', 'normalized', 'position', [0,0,1,1]);
    hold on
    trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'w', 'specularstrength', 0);
    
    if 0%i == 0
        for j = 1:size(e,1)
            l = e(j,:);
            plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'r.-');
        end
    end
    
    if 0%i == 1
        plot3(xyz(v,1), xyz(v,2), xyz(v,3), 'r.');
    end
    
    
    if 0%i == 0
        l = [17 18];
        norm(xyz(l(1),:) - xyz(l(2),:))
        plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'r.-');
        
        x = sum(xyz(l,:),1) / length(l);
        rng = 0.1;
        aabb = x([1,1,2,2,3,3]) + rng*repmat([-1,1],1,3);
    end
    
    axis image vis3d
    %     axis(aabb);
    daspect([1,1,1]);
    %     axis vis3d
    view(3);%-70,-80)
    camproj('persp');
    camlight(30,30);
    
    edg = compute_edges(tri)';
    ledg = sqrt(sum( (xyz(edg(:,1),:) - xyz(edg(:,2),:)).^2, 2));
%     figure;
%     hist(ledg)
    
    fprintf('hmin = %f, hmax = %f, havg = %f\n', min(ledg), max(ledg), mean(ledg));
    
end

return

%%
nsurf = 112;
cls = colorcet( 'I2', 'N', nsurf );
cls = cls(randperm(nsurf),:);
cls = CC( cls, 0.0, 0.8, 1.5 );

figure('units', 'normalized', 'position', [0,0,1,1]);
hold on
for i = [5,10]%1:nsurf
    c = readCoeffs2(sprintf('../cases/jouke/init/coef/c_%3.3d.cheb',i));
    si = surf_chebyshev(c, 1);
    set(si, 'facecolor', cls(i,:) , 'specularstrength', 0);
end

% trisurf(tri, xyz(:,1), xyz(:,2), xyz(:,3), 'facecolor', 'none');
xe = [
    1.0000000000000000      -0.14106846055019359       0.19416407864998739     
   1.0000000000000000      -0.14235402294502350       0.19322352897969691
    ];
xm = [0.99601985123916004        9.1728517453240288E-072   6.9526549333764843E-310];
plot3(xe(:,1), xe(:,2), xe(:,3), 'k.-');
plot3(xm(1), xm(2), xm(3), 'r*');

% for j = 1:size(e,1)
%     l = e(j,:);
%     plot3(xyz(l,1), xyz(l,2), xyz(l,3), 'r.-');
% end

uv = [0.99471843310055408        1.0000000000000000];
c = readCoeffs2(sprintf('../cases/jouke/init/coef/c_010.cheb',i));
ym = ICT2unstr(c, uv);
plot3(ym(1), ym(2), ym(3), 'bo');




axis image vis3d
%     axis(aabb);
daspect([1,1,1]);
%     axis vis3d
view(3);%-70,-80)
camproj('persp');
camlight(30,30);




