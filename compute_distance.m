% Function calculates the coupling strength of Scabrous signal between cell I and cell J
function [cIJ,lD] = compute_distance(latticeX,latticeY,LX,LY,nC)

D = zeros(nC);

[lX,I] = sort(latticeX);

dMax = lX(end) - lX(1);
a = LX;

for k = nC:-1:1
    for l = nC:-1:k+1
        
        dI = abs(latticeX(k) - latticeX(l));
        
        D(k,l) = sqrt((min(dI, abs(a-dI)))^2 + ...
            (latticeY(k) - latticeY(l))^2);
        
    end
end
% This defines cIJ as in equation 4
D = D + D';
lD = 1.75*sqrt((LY*LX)/nC);

D = D./lD;
cIJ = exp(-0.1*(D-3).^2);

cIJ(D == 0 ) = 0;

end