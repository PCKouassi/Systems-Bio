% Function runs the activation wave
function outPut = activation_wave(t,T,y,L,a,l)
    
    if nargin == 4
        a = 1e-2;
        l = 1;
    end
    
    %x = (y./L- t/T);
    %outPut = 0.5*(1-tanh(x/a));
    
    
    outPut = 0.5*(1 - tanh((y - L*(t/T))/(a*l)));
end