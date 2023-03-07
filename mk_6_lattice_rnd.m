% This function makes the lattice of cells

function [nC, vX, vY, c, v, idxNbndr,idxVbndr, latticeX, latticeY] =...
                            mk_6_lattice_rnd(nN,LY,LX,perturbationLevel)


if nargin < 4
    perturbationLevel = 0.2;
end

lN = LX/(nN);
nY = floor(LY/(lN*sin(pi/3)));
nX = nN;

e1 = lN;
e2 = 0.5*lN + 1i*lN*sin(pi/3);



fc = zeros(nX+4,nY+4);

for k = 1:nX+4
    for l = 1:nY+4
        
        fc(k,l) = 0.5*lN*(1+2*1i*tan(pi/6)) + ...
            e1*(k-2 - 0.5*(l-2 + mod(l,2))) + e2*(l-2)+...
                              perturbationLevel*lN*(randn(1) +1i*randn(1));
    end
end


latticeX = real(fc(:));
latticeY = imag(fc(:));


N = length(latticeX);

idxBad = [];
idxOK = [];

idxVok = [];
idxVbad = [];
idxNbndr = [];


[vRand, cRand] = voronoin([latticeX latticeY]);



for n = 1:N
    if (latticeX(n) > 1e-5*LX) && (latticeX(n) <= 1.01*LX) &&...
            (latticeY(n) > 1e-5*LY) && (latticeY(n) <= LY)
        
        idxOK(end+1) = n;    
        idxVok = unique([idxVok,cRand{n}]);
        
        if 1 == 0
            hold on
            patch(vRand([cRand{n}(1:end),cRand{n}(1)],1),...
            vRand([cRand{n}(1:end),cRand{n}(1)],2),...
            [1 0 0])
        end
        
        
    end
end

% number of good cells
nC = length(idxOK);
if nC ~= nX*nY
    fprintf('nC = %d, nX*nY = %d\n',nC, nX*nY)
end

% set of good vertices
vX = vRand(idxVok,1);
vY = vRand(idxVok,2);


% we need to rename veritices to remove unused vertices
vOld2New = zeros(1,size(vRand,1));

for i = 1:length(vX)
    
    vOld2New(idxVok(i)) = i;
end

% vertices that form particular cell
c = cell(1,nC);

for n = 1:nC
    
    c{n} = vOld2New(cRand{idxOK(n)});
end

% boundary vertices
idxVbndr = vOld2New(intersect(idxVok, idxVbad));


for n = 1:nC
    
    if intersect(idxVbndr, c{n})
        idxNbndr(end+1) = n;
    end
end


% list of cells that share each vertix
v = cell(1,length(vX));

for n = 1:nC
    
    for i = c{n}
        v{i}(end+1) = n;
    end
end


latticeX = latticeX(idxOK);
latticeY = latticeY(idxOK);

