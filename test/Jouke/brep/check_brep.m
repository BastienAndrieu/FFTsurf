clc; clear; close all


verts = importdata('verts.dat');
edges = importdata('edges.dat'); % face, orig, prev(2), next(2)
faces = importdata('faces.dat');

nv = size(verts,1);
nh = size(edges,1);
ne = nh/2;
nf = max(edges(:,1));

vxyz = verts(:,1:3);
v2he = round(verts(:,4:5));


cle = hsv(ne);
cle = cle(randperm(ne),:);

clf = colorcet( 'I2', 'N', nf );
clf = clf(randperm(nf),:);

falph = 1;%

%%
figure;
hold on

lv = 1:nv;
% plot3( vxyz(lv,1), vxyz(lv,2), vxyz(lv,3), 'k.');


for iv = 1:nv
    text( vxyz(iv,1), vxyz(iv,2), vxyz(iv,3), num2str(iv), 'color', 'b' );
end

% plot faces
f2v = {};
fb = [];
fracb = 0;%0.05;
i = 1;
for ifa = 1:nf
    f2v{ifa} = [];
    ih = 2*(faces(i,1) - 1) + faces(i,2);
    while 1
        iv = edges(ih,2);
        if ~isempty(f2v{ifa})
            if iv == f2v{ifa}(1); break; end
        end
        f2v{ifa} = [f2v{ifa}, iv];
        next = edges(ih,5:6);
        ih = 2*(next(1) - 1) + next(2);
    end
    i = i + 2 + faces(i+1,1);
    
    b = sum(vxyz(f2v{ifa},:),1) / length(f2v{ifa});
    fb = [fb; b];
    xf = (1 - fracb) * vxyz(f2v{ifa},:) + fracb * repmat(b, length(f2v{ifa}), 1);
    
    
    l = [1:length(f2v{ifa}),1];
    %     plot3( xf(l,1), xf(l,2), xf(l,3), '-', 'color', clf(ifa,:) );
        patch( xf(l,1), xf(l,2), xf(l,3), 'r', ...
            'facecolor', clf(ifa,:), 'edgecolor', 'none', ...0.5*clf(ifa,:), ...
            'facealpha', falph );
end


% plot halfedges as arrows
fracc = 0.08;
fracv = 0.05;
heq = zeros(6,nh);
for ih = 1:nh
    if edges(ih,1) == 0; continue; end
    next = edges(ih,5:6);
    jh = 2*(next(1) - 1) + next(2);
    l = [edges(ih,2), edges(jh,2)];
    
    a = vxyz(l(1),:);
    v = vxyz(l(2),:) - vxyz(l(1),:);
    m = a + 0.5*v;
    
    a = a + fracv*v + fracc * (fb(edges(ih,1),:) - m);
    v = (1 - 2*fracv)*v;
    heq(:,ih) = [a, v]';
    
    quiver3( a(1), a(2), a(3), v(1), v(2), v(3), 0, 'k' );
    
    j = mod(ih,2);
    i = 1 + (ih - j)/2;
    if j == 0
        j = 2;
        i = i - 1;
    end
    text( a(1) + 0.5*v(1), a(2) + 0.5*v(2), a(3) + 0.5*v(3), sprintf('%d,%d',i,j), 'color', 'r');
end

if 0 % plot prev, next & twin of random halfedge
    while 1
        ih = randi(nh);
        if edges(ih,1) ~= 0; break; end
    end
    quiver3( heq(1,ih), heq(2,ih), heq(3,ih), heq(4,ih), heq(5,ih), heq(6,ih), 0, 'k' );
    prev = edges(ih,3:4);
    jh = 2*(prev(1) - 1) + prev(2);
    quiver3( heq(1,jh), heq(2,jh), heq(3,jh), heq(4,jh), heq(5,jh), heq(6,jh), 0, 'r' );
    next = edges(ih,5:6);
    jh = 2*(next(1) - 1) + next(2);
    quiver3( heq(1,jh), heq(2,jh), heq(3,jh), heq(4,jh), heq(5,jh), heq(6,jh), 0, 'g' );
    % twin
    if mod(ih,2) == 0
        jh = ih - 1;
    else
        jh = ih + 1;
    end
    if edges(jh,1) > 0
        quiver3( heq(1,jh), heq(2,jh), heq(3,jh), heq(4,jh), heq(5,jh), heq(6,jh), 0, 'b' );
    end
end

axis image vis3d

axlim = zeros(1,6);
for i = 1:3
    axlim(2*i-1) = min(vxyz(:,i));
    axlim(2*i) = max(vxyz(:,i));
end
rng = axlim(2:2:6) - axlim(1:2:5);
axlim(1:2:5) = axlim(1:2:5) - 0.1*rng;
axlim(2:2:6) = axlim(2:2:6) + 0.1*rng;
axis(axlim); daspect([1,1,1]);

view(3)
camproj('persp');



% lbndedg = [188, 196, 228, 238, 252, 260, 261, 271, 285, 293, 294, 304, 318, 326, 327, 337];
for ie = 1:ne%lbndedg
    for jh = 1:2
        ih = 2*(ie-1) + jh;
        if edges(ih,1) == 0
            ih = 2*(ie-1) + 1 + mod(jh,2);
            next = edges(ih,5:6);
            kh = 2*(next(1) - 1) + next(2);
            l = [edges(ih,2), edges(kh,2)];
            plot3( vxyz(l,1), vxyz(l,2), vxyz(l,3) );
            break
        end
    end
end


% lver = 9;
% plot3( vxyz(lver,1), vxyz(lver,2), vxyz(lver,3), 'k.');
% 
% 
% % ledg = [14,2;
% %     13,1
% %     23,1];
% ledg = [13,1;
%     16,1;
%     15,1;
%     14,1];
% cl = hsv(size(ledg,1));
% for i = 1:size(ledg,1)
%     ih = 2*(ledg(i,1) - 1) + ledg(i,2);
%     quiver3( heq(1,ih), heq(2,ih), heq(3,ih), heq(4,ih), heq(5,ih), heq(6,ih), 0, 'color', cl(i,:) );
% end




%% Feature edges & vertices
lfeat_edge = importdata('feat_edge.dat')';
feat_vert = importdata('feat_vert.dat');

lfeat_vert = find(feat_vert(:,1) == 1)';


%% Hyperfaces & Hyperedges
fid = fopen('hyperfaces.dat','r');
nhf = str2num(fgetl(fid));
for i = 1:nhf
    hyperface(i).nf = str2num(fgetl(fid));
    hyperface(i).faces = str2num(fgetl(fid));
end
fclose(fid);

fid = fopen('hyperedges.dat', 'r');
nhe = str2num(fgetl(fid));
for i = 1:nhe
    hyperedge(i).ne = str2num(fgetl(fid));
    hyperedge(i).verts = str2num(fgetl(fid));
    hyperedge(i).edges = zeros(hyperedge(i).ne,2);
    for j = 1:hyperedge(i).ne
        hyperedge(i).edges(j,:) = str2num(fgetl(fid));
    end
end
fclose(fid);


clhf = colorcet( 'I2', 'N', max(2,nhf) );
% clhf = clhf(randperm(nhf),:);

clhe = 0.6 * colorcet( 'R3', 'N', max(2,nhe) );
% clhe = clhf(randperm(nhe),:);

figure('name','Hyperfaces & Hyperedges');
hold on


if 1 % plot hyperfaces
    for ihf = 1:nhf
        for ifa = hyperface(ihf).faces
            xf = vxyz(f2v{ifa},:);
            l = [1:length(f2v{ifa}),1];
            patch( xf(l,1), xf(l,2), xf(l,3), 'r', ...
                'facecolor', clhf(ihf,:), 'edgecolor', 'none', ...
            'facealpha', falph );
        end
    end
end



% for ih = 1:nh
%     if edges(ih,1) == 0; continue; end
%     j = mod(ih,2);
%     i = 1 + (ih - j)/2;
%     if j == 0
%         j = 2;
%         i = i - 1;
%     end
%     text( ...
%         heq(1,ih) + 0.5*heq(4,ih), ...
%         heq(2,ih) + 0.5*heq(5,ih), ...
%         heq(3,ih) + 0.5*heq(6,ih), ...
%         sprintf('%d,%d',i,j), 'color', 'r');
%     edges(ih,2)
% end




if 1 % plot hyperedges
    for ihe = 1:nhe
        for i = 1:hyperedge(ihe).ne
            iedge = hyperedge(ihe).edges(i,:);
            ih = (2*iedge(1) - 1) + iedge(2);
            next = edges(ih,5:6);
            kh = 2*(next(1) - 1) + next(2);
            l = [edges(ih,2), edges(kh,2)];
%             plot3( vxyz(l,1), vxyz(l,2), vxyz(l,3), '-', 'color', clhe(ihe,:), 'linewidth', 2 );
%             quiver3(vxyz(l(1),1), vxyz(l(1),2), vxyz(l(1),3), ...
%                 vxyz(l(2),1) - vxyz(l(1),1), vxyz(l(2),2) - vxyz(l(1),2), vxyz(l(2),3) - vxyz(l(1),3), ...
%                 0, 'color', clhe(ihe,:), 'linewidth', 1 );
            quiver3( heq(1,ih), heq(2,ih), heq(3,ih), ...
                heq(4,ih), heq(5,ih), heq(6,ih), 0, 'color', clhe(ihe,:), 'linewidth', 1.4 );
        end
    end
else % plot feature edges
    for ie = lfeat_edge
        for jh = 1:2
            ih = 2*(ie-1) + jh;
            if edges(ih,1) ~= 0
                next = edges(ih,5:6);
                kh = 2*(next(1) - 1) + next(2);
                l = [edges(ih,2), edges(kh,2)];
                plot3( vxyz(l,1), vxyz(l,2), vxyz(l,3), 'k-' );
                break
            end
        end
    end
end

plot3( vxyz(lfeat_vert,1), vxyz(lfeat_vert,2), vxyz(lfeat_vert,3), 'k.', 'markersize', 8 );
% for iv = lfeat_vert
%     text(vxyz(iv,1), vxyz(iv,2), vxyz(iv,3), num2str(feat_vert(iv,2)), 'color', 'r');
% end



axis image vis3d
axis(axlim); daspect([1,1,1]);

view(3)
camproj('persp');

return
%%
v2f = {};
fid = fopen('v2f.dat','r');
for ive = 1:nv
    n = str2num(fgetl(fid));
    v2f{ive} = str2num(fgetl(fid));
end
fclose(fid);



figure;
hold on

lv = randi(nv);
plot3( vxyz(lv,1), vxyz(lv,2), vxyz(lv,3), 'k.');

for iv = lv
    ie = v2he(iv,:);
    ih = 2*(ie(1) - 1) + ie(2);
    quiver3( heq(1,ih), heq(2,ih), heq(3,ih), heq(4,ih), heq(5,ih), heq(6,ih), 0, 'k' );
    
    clv2f = colorcet( 'I2', 'N', length(v2f{iv}) );
    %     clv2f = hsv(length(v2f{iv}));
    for jfa = 1:length(v2f{iv})
        ifa = v2f{iv}(jfa);
        xf = (1 - fracb) * vxyz(f2v{ifa},:) + fracb * repmat(fb(ifa,:), length(f2v{ifa}), 1);
        l = [1:length(f2v{ifa}),1];
        patch( xf(l,1), xf(l,2), xf(l,3), 'r', ...
            'facecolor', clv2f(jfa,:), 'edgecolor', 0.5*clv2f(jfa,:), ...
            'facealpha', 0.5 );
        
        %         if jfa < length(v2f{iv})
        %             a = fb(ifa,:);
        %             v = fb(v2f{iv}(jfa+1),:) - fb(ifa,:);
        %             quiver3( a(1), a(2), a(3), v(1), v(2), v(3), 0, 'k');
        %         end
    end
end

axis image vis3d
view(3)
camproj('persp');