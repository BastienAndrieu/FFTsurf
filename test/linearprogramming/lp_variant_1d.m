function [ x, feasible ] = lp_variant_1d( A )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

n = size(A,1);

[ L, R ] = deal( -Inf, +Inf );
bounded = [0,0];

for i = 1:n
    if A(i,1) < 0
        bounded(2) = 1;
        R = min( R, -A(i,2)/A(i,1) );
    elseif A(i,2) > 0
        bounded(1) = 1;
        L = max( L, -A(i,2)/A(i,1) );
    end
end

if L > R
    feasible = 0;
    return
end

feasible = 1;
if all(bounded)
    x = min( abs(L), abs(R) );
elseif bounded(1)
    x = L;
elseif bounded(2)
    x = R;
else
    feasible = -1;
end

end

