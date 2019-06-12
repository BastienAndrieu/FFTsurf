clc; clear; close all

fracc = 0.08;
fracv = 0.05;
fracb = 0;%0.05;
falph = 1;

figure;
for i = 1:2
     subplot(1,2,i)
    hold on
    
    if i == 1
        pth = '../../../init/brep/';
    else
        pth = '';
    end
    
    brep(i).verts = importdata([pth, 'verts.dat']);
    brep(i).edges = importdata([pth, 'edges.dat']);
    brep(i).faces = importdata([pth, 'faces.dat']);
    
    nf = max(brep(i).edges(:,1));
    
    clf = colorcet('I2', 'N', nf);
    
    vxyz = brep(i).verts(:,1:3);
    
    fb = [];
    f2v = {};
    j = 1;
    for ifa = 1:nf
        f2v{ifa} = [];
        ih = 2*(brep(i).faces(j,1) - 1) + brep(i).faces(j,2);
        while 1
            iv = brep(i).edges(ih,2);
            if ~isempty(f2v{ifa})
                if iv == f2v{ifa}(1); break; end
            end
            f2v{ifa} = [f2v{ifa}, iv];
            next = brep(i).edges(ih,5:6);
            ih = 2*(next(1) - 1) + next(2);
        end
        j = j + 2 + brep(i).faces(j+1,1);
        
        b = sum(vxyz(f2v{ifa},:),1) / length(f2v{ifa});
        fb = [fb; b];
        
        xf = (1 - fracb) * vxyz(f2v{ifa},:) + fracb * repmat(b, length(f2v{ifa}), 1);
        l = [1:length(f2v{ifa}),1];
        
%         patch( xf(l,1), xf(l,2), xf(l,3), 'r', ...
%             'facecolor', clf(ifa,:), 'edgecolor', 'none', ...0.5*clf(ifa,:), ...
%             'facealpha', falph );
    end
    
    
    
  
    brep(i).heq = zeros(6,size(brep(i).edges,1));
    for ih = 1:size(brep(i).edges,1)
        if brep(i).edges(ih,1) == 0; continue; end
        next = brep(i).edges(ih,5:6);
        jh = 2*(next(1) - 1) + next(2);
        l = [brep(i).edges(ih,2), brep(i).edges(jh,2)];
        
        a = vxyz(l(1),:);
        v = vxyz(l(2),:) - vxyz(l(1),:);
        m = a + 0.5*v;
        
        a = a + fracv*v + fracc * (fb(brep(i).edges(ih,1),:) - m);
        v = (1 - 2*fracv)*v;
        heq(:,ih) = [a, v]';
    end
    
    
    fid = fopen([pth, 'hyperfaces.dat'],'r');
    nhf = str2num(fgetl(fid));
    for j = 1:nhf
        hypg(i).hyperface(j).nf = str2num(fgetl(fid));
        hypg(i).hyperface(j).faces = str2num(fgetl(fid));
    end
    fclose(fid);

    fid = fopen([pth, 'hyperedges.dat'], 'r');
    nhe = str2num(fgetl(fid));
    for j = 1:nhe
        hypg(i).hyperedge(j).ne = str2num(fgetl(fid));
        hypg(i).hyperedge(j).verts = str2num(fgetl(fid));
        hypg(i).hyperedge(j).hyperfaces = str2num(fgetl(fid));
        hypg(i).hyperedge(j).edges = zeros(hypg(i).hyperedge(j).ne,2);
        for k = 1:hypg(i).hyperedge(j).ne
            hypg(i).hyperedge(j).edges(k,:) = str2num(fgetl(fid));
        end
    end
    fclose(fid);
    
    clhe = 0.6 * colorcet( 'R3', 'N', max(2,nhe) );
    for ihe = 1:nhe
        for j = 1:hypg(i).hyperedge(ihe).ne
            iedge = hypg(i).hyperedge(ihe).edges(j,:);
            ih = 2*(iedge(1) - 1) + iedge(2);
            next = brep(i).edges(ih,5:6);
            kh = 2*(next(1) - 1) + next(2);
            l = [brep(i).edges(ih,2), brep(i).edges(kh,2)];
            if j == 1
                plot3(heq(1,ih), heq(2,ih), heq(3,ih), '*', 'color', clhe(ihe,:));
                 a = heq(1:3,ih) + 0.5*heq(4:6,ih);
            text(a(1), a(2), a(3), num2str(ihe), 'color', clhe(ihe,:), 'backgroundcolor', 'w');
            end
            quiver3( heq(1,ih), heq(2,ih), heq(3,ih), ...
                heq(4,ih), heq(5,ih), heq(6,ih), 0, ...
                'color', clhe(ihe,:), 'linewidth', 1.4, 'MaxHeadSize', 10 );
        end
    end
    
    
%     plot3(brep(i).verts(:,1), brep(i).verts(:,2), brep(i).verts(:,3), '*');
    
    for j = 1:size(brep(i).verts,1)
        x = brep(i).verts(j,:);
        text(x(1), x(2), x(3), num2str(j), 'backgroundcolor', 'w');
    end
    
    axis image
    view(3)
end


