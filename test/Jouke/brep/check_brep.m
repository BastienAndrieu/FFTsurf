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


%%
figure;
hold on

lv = 1:nv;
plot3( vxyz(lv,1), vxyz(lv,2), vxyz(lv,3), 'k.');


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
        patch( xf(l,1), xf(l,2), xf(l,3), 'r', 'facecolor', clf(ifa,:), 'edgecolor', 0.5*clf(ifa,:) );
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
end

if 0
    % plot prev, next & twin of random halfedge
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
        patch( xf(l,1), xf(l,2), xf(l,3), 'r', 'facecolor', clv2f(jfa,:), 'edgecolor', 0.5*clv2f(jfa,:) );
        
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