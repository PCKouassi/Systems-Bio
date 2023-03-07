% Runs equation 2
function aUU = Du(u,a0)

if nargin == 1
    
    a0 = 0.05;
end

aUU = u.*(a0 + (1 - a0).*(2*u.^3)./(1 + u.^2));
%aUU = 0 ;