clc; clear; close all
falph = 1;%0.5;%


folders = {'offset_3'; 'offset_12'};

figure('units', 'normalized', 'position', [.15, .15, .7, .7]);

for i = 1:2
    verts = importdata([folders{i},'/verts.dat']);
    instant(i).edges = importdata([folders{i},'/edges.dat']);
    instant(i).faces = importdata([folders{i},'/faces.dat']);
    
    instant(i).nv = size(verts,1);
    nh = size(instant(i).edges,1);
    instant(i).ne = nh/2;
    instant(i).nf = max(instant(i).edges(:,1));
    
    instant(i).vxyz = verts(:,1:3);
    instant(i).v2he = round(verts(:,4:5));
    
    instant(i).lfeat_edge = importdata([folders{i},'/feat_edge.dat'])';
    feat_vert = importdata([folders{i},'/feat_vert.dat']);
    instant(i).lfeat_vert = find(feat_vert(:,1) == 1)';
    
    fid = fopen([folders{i},'/hyperfaces.dat'],'r');
    instant(i).nhf = str2num(fgetl(fid));
    for j = 1:instant(i).nhf
        instant(i).hyperface(j).nf = str2num(fgetl(fid));
        instant(i).hyperface(j).faces = str2num(fgetl(fid));
    end
    fclose(fid);
    
    instant(i).f2hf = zeros(1,instant(i).nf);
    for j = 1:instant(i).nhf
        for k = instant(i).hyperface(j).faces
            instant(i).f2hf(k) = j;
        end
    end
    
    fid = fopen([folders{i},'/hyperedges.dat'], 'r');
    instant(i).nhe = str2num(fgetl(fid));
    for j = 1:instant(i).nhe
        instant(i).hyperedge(j).ne = str2num(fgetl(fid));
        instant(i).hyperedge(j).verts = str2num(fgetl(fid));
        instant(i).hyperedge(j).hyperfaces = str2num(fgetl(fid));
        for k = 1:instant(i).hyperedge(j).ne
            instant(i).hyperedge(j).edges(k,:) = str2num(fgetl(fid));
        end
    end
    fclose(fid);
    
%     instant(i).he2hf = zeros(instant(i).nhe,2);
%     for j = 1:instant(i).nhe
%         for k = 1:instant(i).hyperedge(j).ne
%                 ihedg = instant(i).hyperedge(j).edges(k,:);
%                 for h = 1:2
%                     ih = 2*(ihedg(1) - 1) + ihedg(2);
%                     f = instant(i).f2hf(instant(i).edges(ih,1));
%                     if instant(i).he2hf(j,h) == 0
%                         instant(i).he2hf(j,h) = f;
%                     end
%                     ihedg(2) = 1 + mod(ihedg(2),2);
%                 end
%         end
%     end
    
    
    instant(i).f2v = {};
    l = 1;
    for ifa = 1:instant(i).nf
        instant(i).f2v{ifa} = [];
        ih = 2*(instant(i).faces(l,1) - 1) + instant(i).faces(l,2);
        while 1
            iv = instant(i).edges(ih,2);
            if ~isempty(instant(i).f2v{ifa})
                if iv == instant(i).f2v{ifa}(1); break; end
            end
            instant(i).f2v{ifa} = [instant(i).f2v{ifa}, iv];
            next = instant(i).edges(ih,5:6);
            ih = 2*(next(1) - 1) + next(2);
        end
        l = l + 2 + instant(i).faces(l+1,1);
    end
    
    
    clhf = colorcet( 'I2', 'N', max(2,instant(i).nhf) );
    clhe = 0.6 * colorcet( 'R3', 'N', max(2,instant(i).nhe) );
    clhe = clhe(randperm(size(clhe,1)),:);
    
    subplot(1,2,i);
    hold on
    
    if 1 % plot hyperfaces
        for ihf = 1:instant(i).nhf
            for ifa = instant(i).hyperface(ihf).faces
                xf = instant(i).vxyz(instant(i).f2v{ifa},:);
                l = [1:length(instant(i).f2v{ifa}),1];
                patch( xf(l,1), xf(l,2), xf(l,3), 'r', ...
                    'facecolor', clhf(ihf,:), 'edgecolor', 0.8*clhf(ihf,:),...'none', ...
                    'facealpha', falph );
            end
        end
    end
    
    if 1 % plot hyperedges
        for ihe = 1:instant(i).nhe
            for j = 1:instant(i).hyperedge(ihe).ne
                iedge = instant(i).hyperedge(ihe).edges(j,:);
                ih = 2*(iedge(1) - 1) + iedge(2);
                next = instant(i).edges(ih,5:6);
                kh = 2*(next(1) - 1) + next(2);
                l = [instant(i).edges(ih,2), instant(i).edges(kh,2)];
                plot3( instant(i).vxyz(l,1), instant(i).vxyz(l,2), instant(i).vxyz(l,3), ...
                    '-', 'color', clhe(ihe,:), 'linewidth', 2 );
            end
        end
    end
    
    if 1 % plot feature vertices
        plot3( ...
            instant(i).vxyz(instant(i).lfeat_vert,1), ...
            instant(i).vxyz(instant(i).lfeat_vert,2), ...
            instant(i).vxyz(instant(i).lfeat_vert,3), 'k.', 'markersize', 10 );
    end
    
%     axis image vis3d
axis([-0.1, 1.25, 0, 0.5, -0.5, 0.5])
daspect([1,1,1])
    axis vis3d
    view(3)
    camproj('persp');
    
    title( folders{i}, 'interpreter', 'none' );
    
    
    fprintf('----------- Instant %d -----------\n',i);
    fprintf('%d verts\n',instant(i).nv)
    fprintf('%d edges\n',instant(i).ne)
    fprintf('%d faces\n',instant(i).nf);
    fprintf('%d hyperedges\n',instant(i).nhe);
    fprintf('%d hyperfaces\n',instant(i).nhf);
    fprintf('\n\n');
    fprintf('hyperface -> faces :\n');
    for ihf = 1:instant(i).nhf
        fprintf('%d :\t',ihf);
        for j = instant(i).hyperface(ihf).faces
            fprintf('%d ',j);
        end
        fprintf('\n');
    end
    fprintf('\n');
    fprintf('hyperedge -> hyperfaces :\n');
    for ihe = 1:instant(i).nhe
        fprintf('%d :\t',ihe);
        for j = 1:2
            fprintf('%d ', instant(i).hyperedge(ihe).hyperfaces(j));
        end
        fprintf('\n');
    end
    fprintf('----------------------------------\n\n\n');
    
    
    
end

if instant(1).nf == instant(2).nf
    f_old2new = 1:instant(1).nf;
else
    f_old2new = zeros(1,instant(1).nf);
end


hf_old2new = zeros(1,instant(1).nhf);
for i = 1:instant(1).nhf
    for j = 1:instant(2).nhf
        lenij = length(intersect(...
            f_old2new(instant(1).hyperface(i).faces), ...
            instant(2).hyperface(j).faces));
        
        if lenij > 0
            hf_old2new(i) = j;
            break
        end
    end
end

hf_old2new

he_old2new = zeros(1,instant(1).nhe);
for i = 1:instant(1).nhe
    hfi = hf_old2new(instant(1).hyperedge(i).hyperfaces);
    for j = 1:instant(2).nhe
        hfj = instant(2).hyperedge(j).hyperfaces;
        if hfi(1) == hfj(1) && hfi(2) == hfj(2)
            he_old2new(i) = j;
        elseif hfi(1) == hfj(2) && hfi(2) == hfj(1)
            he_old2new(i) = -j;
        end
    end
end

he_old2new

