clc; clear; close all

n = 7;
x = linspace(0,1,n)';
y = rand(n,1);

b = [x,y];

fid = fopen( 'b.dat', 'w' );
fprintf( fid, '%d %d\n', size(b,1), size(b,2) );
for j = 1:size(b,2)
    for i = 1:size(b,1)
        fprintf( fid, '%16.15e\n', b(i,j) );
    end
end
fclose(fid);

t = rand
fid = fopen( 't.dat', 'w' );
fprintf( fid, '%16.15e\n', t );
fclose(fid);