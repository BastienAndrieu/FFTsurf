clc; clear; close all


[ a14, a24, a34 ] = deal( 5 );
[ ai1, ai2, ai3, ai4 ] = deal( 0.3, -0.7, -1.0, 4.0 );
[ aj1, aj2, aj3, aj4 ] = deal( -2.0, -1.0, 1.0, 6.0 );
[ ak1, ak2, ak3, ak4 ] = deal( 0.5, -0.6, 1.0, -3.0 );
[ al1, al2, al3, al4 ] = deal( -0.1, 0.2, -1.0, 3.0 );
A = [
    [ -eye(3) , [a14 ; a24; a34] ];
    ai1, ai2, ai3, ai4 ;
    aj1, aj2, aj3, aj4 ;
    ak1, ak2, ak3, ak4 ;
    al1, al2, al3, al4 ;
    ];
A(4:end,:) = A(3+randperm(size(A,1)-3),:);
Ap = [ eye(3), zeros(3,1) ];

% c = [-1.0;-1.0;-1.0];
% c = -rand(3,1);
c = [
    -0.999080394761361
    -0.171121066356432
    -0.032600820530528];
% c = [1.0;-1.0;0.0];
% c = [-1.0;0.0;1.0];

x = -A(1:3,1:3) \ A(1:3,4);



% Ai = [
%     [ -eye(2) , [a14 ; a24] ];
%     ai1/ai3, ai2/ai3, a34 + ai4/ai3 ];
% Api = [ [eye(2), zeros(2,1)];
%      -ai1/ai3, -ai2/ai3, -ai4/ai3 ];
% 
% ci = [ c(1) - ai1/ai3*c(3) ; c(2) - ai2/ai3*c(3) ];
% visu_lpp( [Ai;Api], ci )
% keyboard

A = [ A(1:3,:); Ap; A(4:end,:) ];
[ x, feasible ] = lpp_solve( A, c, x )






% Av = [A ; Ap];
Av = A;
v = con2vert( -Av(:,1:3), Av(:,4) );
fv = v*c;
q = convhulln( v );

figure;
hold on
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
colormap(colorcet( 'L18', 'N', 16 ));

if feasible
    plot3( x(1), x(2), x(3), 'b*');
end


view(145,30)
axis( [0, a14, 0, a24, 0, a34] ); % axis image
axis vis3d
daspect([1,1,1])
grid on

xlabel('x');
ylabel('y');
zlabel('z');