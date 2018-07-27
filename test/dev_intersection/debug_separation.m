clc; clear; close all

addpath('/stck/bandrieu/Bureau/CYPRES/Intersections/separation/');
addpath('/stck/bandrieu/Bureau/CYPRES/FFTsurf/FORTRAN/Intersection/Curve-Surface_2/');

x = [
    4.7846698502724474E-017   1.4870700157376682E-014   1.6028844918025698E-015
%     -1.2442823401414761E-016  -4.4524056019758923E-009   1.4953316362920077E-015
%     -1.2442823388657825E-016   1.9394666090147283E-009   1.5005358067199381E-015
%     1.6423025911719281E-008  -8.2718061255302767E-025  -5.2041706347255744E-017
    ];

for i = 1:2
    s(i).xyz = importdata(sprintf('sep%d.dat',i));
    if i == 1
        s(i).b = readCoeffs1(sprintf('sep%d_poly.bern',i));
        sqrt(sum((s(i).b([1,end],:) - repmat(x,2,1)).^2,2))
    else
        s(i).b = readCoeffs2(sprintf('sep%d_poly.bern',i));
        q = [
            reshape(s(i).b(1,1,:)    , 1, 3);
            reshape(s(i).b(end,1,:)  , 1, 3);
            reshape(s(i).b(1,end,:)  , 1, 3);
            reshape(s(i).b(end,end,:), 1, 3)
            ];
        sqrt(sum((q - repmat(x,4,1)).^2,2))
    end
end

cl = ['r','b'];

for j = 2%1:2
    figure;
    hold on
    
    for i = 1:2
        if j == 1
            plot3(s(i).xyz(:,1), s(i).xyz(:,2), s(i).xyz(:,3), '.', 'color', cl(i));
        else
            if size(s(i).b,3) == 3
                surf(s(i).b(:,:,1) - x(1), s(i).b(:,:,2) - x(2), s(i).b(:,:,3) - x(3), ...
                    'facecolor', 'none', 'edgecolor', cl(i));
            else
                plot3(s(i).b(:,1) - x(1), s(i).b(:,2) - x(2), s(i).b(:,3) - x(3), ...
                    '.-', 'color', cl(i));
            end
        end
    end
    
    if j == 1
        plot3(0, 0, 0, 'kx');
        [ vec, stat ] = separating_plane( s(1).xyz, s(2).xyz )
        
%         s(1).xyz*vec, s(2).xyz*vec
        
        if stat == 0
            plot_separating_plane(vec);
        end
        axis( repmat( [-1,1], 1, 2 ) );
        daspect([1,1,1]);
    else
        plot3(x(1), x(2), x(3), 'k*');
        axis image;
    end
%     axis vis3d
%     view(3)
    view(120,30)
    camlight(30,30);
end