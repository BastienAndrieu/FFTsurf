clc; clear; %close all

instant = 115;

tri = importdata('../connect.dat')';
xyz = importdata(sprintf('../pos_%3.3d.dat',instant))';

nt = size(tri,2);
nv = size(xyz,2);

norv = zeros(3,nv);
v2t = {};
for iv = 1:nv
    v2t{iv} = [];
end

% compute pseufor-normals at vertices
for it = 1:nt
    vec = cross( ...
        xyz(:,tri(2,it)) - xyz(:,tri(1,it)), ...
        xyz(:,tri(3,it)) - xyz(:,tri(1,it)) );
    for jv = 1:3
        iv = tri(jv,it);
        norv(:,iv) = norv(:,iv) + vec;
        v2t{iv} = [v2t{iv}, it];
    end
end

% compute discrete minimum curvature radius at vertices
rad = 1e6 * ones(nv,1);
v2v = {};
for iv = 1:nv
    norv(:,iv) = norv(:,iv) / norm(norv(:,iv));
    
    lv = tri(:,v2t{iv});
    lv = setdiff(unique(lv(:)'), iv);
    v2v{iv} = lv;
    
    for jv = v2v{iv}
        vec = xyz(:,jv) - xyz(:,iv);
        rad(iv) = min(rad(iv), sum(vec.^2)/abs(dot(vec, norv(:,iv))));
    end
end


ratio = 1.0/0.33;
ec = 'none';%'k';%

hve = rad / min(rad);
hve = min(ratio, hve);
% if max(hve) > min(hve) + eps('double')
%     hve = (hve - min(hve)) / (max(hve) - min(hve));
% end

% figure('units','normalized', 'position', [.1 .11 .7 .7]);
% trisurf(tri', xyz(1,:), xyz(2,:), xyz(3,:), hve, ...
%     'facecolor', 'interp', 'edgecolor', ec);
% axis image vis3d
% view(3);
% camproj('persp');
% colormap(flipud(cbrewer('seq','OrRd',256)))


% smoothing
npass = 30;
tmp = zeros(nv,1);
mdrx = [];
sumf = [];
for ipass = 1:npass
    max_dr_dx = 0;
    tmp(1:nv) = 0.0;
    for ivert = 1:nv
        sumw = 0.0;
        for jvert = v2v{ivert}
            w = 1.0 / norm(xyz(:,jvert) - xyz(:,ivert));
            maxdr = max(max_dr_dx, abs(hve(jvert) - hve(ivert))*w);
            sumw = sumw + w;
            tmp(ivert) = tmp(ivert) + w * hve(jvert);
        end
        w = sumw / length(v2v{ivert});
        tmp(ivert) = tmp(ivert) + w * hve(ivert);
        sumw = sumw + w;
        tmp(ivert) = tmp(ivert) / sumw;
    end
    hve(1:nv) = tmp(1:nv);
%     if ipass > 1
%         hve(1:nv) = hve(1:nv) * sumf(1)/sum(hve(1:nv));
%     end
    mdrx = [mdrx; maxdr];
    sumf = [sumf; sum(hve(1:nv))];
end


figure('name','h/h_{min}','units','normalized', 'position', [.1 .11 .7 .7]);
trisurf(tri', xyz(1,:), xyz(2,:), xyz(3,:), hve, ...
    'facecolor', 'interp', 'edgecolor', ec);
axis image vis3d
view(3);
camproj('persp');
colormap(flipud(cbrewer('seq','Blues',256)))
colorbar

% figure;
% plot(mdrx,'.-');

% figure;
% plot(sumf/sumf(1),'.-');




area_act = zeros(nt,1);
area_ide = zeros(nt,1);
for i = 1:nt
    u = xyz(:,tri(2,i)) - xyz(:,tri(1,i));
    v = xyz(:,tri(3,i)) - xyz(:,tri(1,i));
    area_act(i) = norm(cross(u,v));
    
    htri = sum(hve(tri(:,i)));
    area_ide(i) = htri^2;
end
area_ide = area_ide * sum(area_act)/sum(area_ide);


figure('name','A_{act}/A_{ide}','units','normalized', 'position', [.1 .11 .7 .7]);
trisurf(tri', xyz(1,:), xyz(2,:), xyz(3,:), area_act./area_ide, 'edgecolor', ec);
axis image vis3d
view(3);
camproj('persp');
colormap(cbrewer('seq','OrRd',32))
colorbar
