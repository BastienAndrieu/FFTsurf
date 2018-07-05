function [  ] = visu_lpp( A, c, x )
%VISU_LPP Summary of this function goes here
%   Detailed explanation goes here


n = size(A,1)
dim = length(c)

if dim > 3
    error('dim > 3 not supported')
end


figure;
hold on

if dim == 1
%     L = -Inf;
%     R = +Inf;
%     
%     np = 0;
%     bounded_LR = [0,0];
%     
%     m = 0; j = 0;
%     for i = 1:n
%         if A(i,1) == 0; continue; end
%         j = j + 1;
%         p = -A(i,2)/A(i,1);
%         m = m + p;
%         if A(i,1) < 0
%             plot( p, 0, 'k>' );
%             np = np + 1;
%             R = min(R, p);
%             bounded_LR(2) = 1;
%         elseif A(i,1) > 0
%             plot( p, 0, 'k<' );
%             L = max(L, p);
%             bounded_LR(1) = 1;
%         end
%     end
%     m = m / j;
%     quiver( m, 1, sign(c), 0 );
p = [];
for i = 1:n
    if A(i,1) == 0; continue; end
    p = [ p , -A(i,2)/A(i,1) ];
end
plim = [ min(p), max(p) ];
for i = 1:n
    if A(i,1) == 0; continue; end
    plot( plim, A(i,1)*plim + A(i,2) );
end

quiver( 0.5*sum(plim), 0, sign(c), 0 );
if nargin > 2
    plot( x, 0, 'b*');
end
    
else

    v = con2vert( -A(:,1:dim), A(:,dim+1) );
    fv = v*c;
%     if nargin > 2
%         if dot(x,c) > min( fv ) + eps('single');
%             min(fv), dot(x,c)
%             fprintf('x is not an optimal point! :(\n')
%             [~,imin] = min( fv );
%             x, v(imin,:)
%             if dim == 2
%                 plot( v(:,1), v(:,2), 'k.' );
%                 plot( v(imin,1), v(imin,2), 'ko' );
%             elseif dim == 3
%                 plot3( v(:,1), v(:,2), v(:,3), 'k.' );
%                 plot3( v(imin,1), v(imin,2), v(imin,3), 'ko' );
%             end
%         end
%     end
    q = convhulln( v );
    
    if dim == 2
        
        if ~isempty(v)
            h = delaunay(v);
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
        
        if nargin > 2
            plot( x(1), x(2), 'b*');
        end
        view(2)
        
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
        
        if nargin > 2
            plot3( x(1), x(2), x(3), 'b*');
        end
        view(145,30)
    end
    
%     view(dim)
end

colormap(colorcet( 'L18', 'N', 16 ));

axis image vis3d
grid on





end

