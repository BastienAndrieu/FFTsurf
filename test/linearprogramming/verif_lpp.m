clc; clear; close all
addpath('/stck/bandrieu/Bureau/CYPRES/LinearProgramming/');


fid = fopen( 'result.dat', 'r' );
fname = strtrim( fgetl(fid) );
stat = str2num( fgetl(fid) );
if stat == 0
    x = str2num( fgetl(fid) );
end
fclose(fid);


fid = fopen( fname, 'r' );
fgetl(fid); %'dimension'
dim = str2num( fgetl(fid) );
fgetl(fid); %'objectif'
c = str2num( fgetl(fid) )';
fgetl(fid); %'contraintes'
nc = str2num( fgetl(fid) );
A = zeros(nc,dim+1);
for i = 1:nc
    A(i,:) = str2num( fgetl(fid) );
end
fclose(fid);


%% Visualisation
% [ sol, feasible ] = solve_bounded_LPP_d( A, c );
% order = [1:3,3+randperm(nc-3)];
% order = [1:6,6+randperm(nc-6)];

% A = A(order,:);
sol = -A(1:dim,1:dim) \ A(1:dim,dim+1);
[ sol, feasible ] = lpp_solve( A, c, sol );


A = [ A ; [ eye(dim), 200*ones(dim,1) ] ];

v = con2vert( -A(:,1:dim), A(:,dim+1) );

if ~isempty(v)
    fv = v*c;
    q = convhulln( v );
    
    figure;
    hold on
    
    if dim == 2
        
        
        if ~isempty(v)
            h = delaunayn(v);
            patch('Faces',h, 'Vertices', v, ...
                'FaceColor','interp',...
                'FaceVertexCData',fv,...
                'CDataMapping','scaled',...
                'EdgeColor','none');
            %     colormap(cbrewer('seq','Blues',16));
            
            for j = 1:size(q,1)
                plot( v(q(j,:),1), v(q(j,:),2), 'k-' );
            end
        end
        
        if stat == 0
            plot( x(1), x(2), 'b*');
        end
        
        if feasible
            plot( sol(1), sol(2), 'bo');
        end
        
        
    elseif dim == 3
        
        patch('Faces',q, 'Vertices', v, ...
            'FaceColor','interp',...
            'FaceVertexCData',fv,...
            'CDataMapping','scaled',...
            'EdgeColor','none');
        %
        [ e2v, e2t ] = triangles2edges( q );
        nor = zeros( size(q,1), 3 );
        for i = 1:size(q,1)
            nor(i,:) = cross( v(q(i,2),:) - v(q(i,1),:), v(q(i,3),:) - v(q(i,1),:) );
            nor(i,:) = nor(i,:) / norm( nor(i,:) );
        end
        for i = 1:size(e2t,1)
            if e2t(i,2) > 0
                if dot( nor(e2t(i,1),:),nor(e2t(i,2),:) ) > 1.0 - 1e-6
                    continue
                end
            end
            plot3( v(e2v(i,:),1), v(e2v(i,:),2), v(e2v(i,:),3), 'k-' );
        end
        %
        if stat == 0
            plot3( x(1), x(2), x(3), 'b*');
        end
        
        if feasible
            plot3( sol(1), sol(2), sol(3), 'bo');
        end
        
    end
    
    colormap(colorcet( 'L18', 'N', 32 ));
    
    view(dim)
    axis image vis3d
    grid on
end