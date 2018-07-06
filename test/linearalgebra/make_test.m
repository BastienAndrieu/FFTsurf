clc; clear; close all

[m,n] = deal(100);

A = 2*rand(m,n) - 1;


fid = fopen('A.dat','w');
fprintf( fid, '%d %d\n', m, n );
for j = 1:n
    for i = 1:m
        fprintf( fid, '%16.15e\n', A(i,j) );
    end
end
fclose(fid);


[U,W,V] = svd(A);

cond(A)



b = rand(m,1);
fid = fopen('b.dat','w');
for i = 1:m
    fprintf( fid, '%16.15e\n', b(i) );
end
fclose(fid);


x = A \ b;
max( abs( A*x - b ) )
