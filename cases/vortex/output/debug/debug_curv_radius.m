clc; clear; close all

instant = 115;


nfaces = 6;
ncg = 100;
ucg = linspace(-1,1,ncg)';


minrcm = Inf;
for ifa = 1:nfaces
    if instant == 0
        c = readCoeffs2(sprintf('../../init/coef/c_%3.3d.cheb',ifa));
    else
        c = readCoeffs2(sprintf('instant_%3.3dc_%3.3d.cheb',instant,ifa));
    end
    
    face(ifa).c = c;
    
    D = cheb_diff2(c);
    D2u = cheb_diff2(D(:,:,:,1));
    D2v = cheb_diff2(D(:,:,:,2));
    
    s = ICT2tpgrid(c,ucg,ucg);
    su = ICT2tpgrid(D(:,:,:,1),ucg,ucg);
    sv = ICT2tpgrid(D(:,:,:,2),ucg,ucg);
    suu = ICT2tpgrid(D2u(:,:,:,1),ucg,ucg);
    suv = ICT2tpgrid(D2u(:,:,:,2),ucg,ucg);
    svv = ICT2tpgrid(D2v(:,:,:,2),ucg,ucg);
    
    cE = sum(su.*su,3);
    cF = sum(su.*sv,3);
    cG = sum(sv.*sv,3);
    
    cW2 = cE.*cG - cF.*cF;
    cW = sqrt(cW2);
    
    sn = zeros(ncg,ncg,3);
    for i = 1:3
        j = 1+mod(i,3);
        k = 1+mod(j,3);
        
        sn(:,:,i) = (su(:,:,j).*sv(:,:,k) - su(:,:,k).*sv(:,:,j))./cW;
    end
    
    cL = sum(sn.*suu,3);
    cM = sum(sn.*suv,3);
    cN = sum(sn.*svv,3);
    
    cK = (cL.*cN - cM.*cM)./cW2;
    cH = 0.5*(cE.*cN + cG.*cL - 2.0*cF.*cM)./cW2;
    
    cS = sqrt(cH.*cH - cK);
    k1 = cH + cS;
    k2 = cH - cS;
    
    rcm = 1.0./max(abs(k1),abs(k2));
    minrcm = min(minrcm, min(rcm(:)));
    
    face(ifa).rcm = rcm;
end

ratio = 1.0/0.33;

maxhve = 0;
minhve = Inf;
for ifa = 1:nfaces
    rcm = face(ifa).rcm / minrcm;
    rcm = min(ratio, rcm);

    maxhve = max(maxhve, max(rcm(:)));
    minhve = min(minhve, min(rcm(:)));
    face(ifa).rcm = rcm;
end

% for ifa = 1:nfaces
%     face(ifa).rcm = (face(ifa).rcm - minhve) / (maxhve - minhve);
% end


figure('units','normalized', 'position', [.1 .11 .7 .7]);
hold on
for i = 1:nfaces
    s = ICT2tpgrid(face(i).c,ucg,ucg);
    surf(s(:,:,1), s(:,:,2), s(:,:,3), face(i).rcm, ...
        'facecolor', 'interp', 'edgecolor', 'none');
end

axis image vis3d
view(3);
camproj('persp');

colormap(flipud(cbrewer('seq','Blues',256)))
colorbar
