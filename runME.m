%% Parameter

% Number of cells in a row
nN = 2^3;

% Size of the unit rectangle of cells
LX = 20; %dont change this

LY = 3*LX;

% Noise level of u. You can play with this parameter
alphaNoise = 1;

% Defines how disordered the cells pack together. zero is perfectly ordered 
latticePerturbationLevel = .2;

% Steepness of the activation wave front
aL = 0.1;

% Speed of the activation wave 
aT = 1;

% Duration of the simulation
T = 10;

% The rate of the activation loops in the ODE. don't change
tauInv = 2;



%%  Cell lattice intiation
[nC, vX, vY, c, v, idxNbndr,idxVbndr,latticeX,latticeY] = ...
    mk_6_lattice_rnd(nN,LY,LX,latticePerturbationLevel);

% Random assignment of state u to each cell I
uI = mod(abs(randn(nC,1)),1);

%% Plot of the initial cell lattice

close all
h = figure(1);


set(gcf, 'Units', 'Normalized', 'OuterPosition',...
    [0.1, 0.04, min(LX/LY,0.9), min(0.96,LY/LX)]);

for nI = 1:nC
    
    patch(  vX([c{nI}(1:end),c{nI}(1)])    ,...
        vY([c{nI}(1:end),c{nI}(1)])    ,...
        [1 uI(nI)^0.3 0]                );
end

drawnow;
axis([0 LX 0 LY])
axis equal
colormap autumn
colorbar


%% Matrix of cell-to-cell distances

[cIJ,lD] = compute_distance(latticeX,latticeY,LX,LY,nC);

%% Set and run the model
tic;

% Set frame size
h = figure(1);
set(gcf, 'Units', 'Normalized', 'OuterPosition',...
    [0.1, 0.04, min(LX/LY,0.9), min(0.96,LY/LX)]);
axis([0 LX 0 LY])
axis equal
setappdata(gca,'LegendColorbarManualSpace',1);
setappdata(gca,'LegendColorbarReclaimSpace',1);


% Set video frame rate
videoFps = 6;
framePerImage = 1;

% Give name for video
writerObj = VideoWriter(sprintf(...
    'output/%02d_fps_rand_%02d_cell_T_%1.1f_a_%02.2f_noise_%02.2f.mp4',...
    videoFps,nN,aT,aL,alphaNoise),'MPEG-4');
writerObj.FrameRate = videoFps;

open(writerObj);

clear('Frame')

% Initial randomly assign values of state u for each cell I
uI = mod(abs(randn(nC,1)),1);

k = 1;
tStart = toc;

% Assign size of discretization step in simulation. this parameter sshould be relatively
%small
dt = 0.1;
 % Define total time simulation runs. The number 2 means it runs twice as long as it takes the wave to pass across the entire field. 
  % This number can be varied but should be greater than 1. 
for tI = 0:dt:max(2*aT*T,T)
 % Define s as in equation 3   
    sI = cIJ*Du(uI);
 % Define f as in equations 6 and 7   
    fI = sigmoidTh(2.*(uI-sI));
 % Define u as in equation 5   
    uI = uI +...
        tauInv*dt*(activation_wave(tI,aT*T,latticeY,LY,aL,lD).*fI - uI)...
        + tauInv*sqrt(dt)*alphaNoise*randn(nC,1)*sqrt(2e-5);
    
    uI = min(max(uI,0),1);
    
    fprintf('t = %2.1f\n',tI)
    
     % nI is the index of cells 
    for nI = 1:nC
        
        patch(vX([c{nI}(1:end),c{nI}(1)]), vY([c{nI}(1:end),c{nI}(1)]),...
            [1 uI(nI)^0.2 0])
    end
     % Assigns color to the value of u in each cell. Yellow is maximal and
     % red is minimal
    drawnow;
    set(gcf, 'Units', 'Normalized', 'OuterPosition',...
        [0.1, 0.04, min(LX/LY,0.9), min(0.96,LY/LX)]);
    axis([0 LX 0 LY])
    axis equal
    
    Frame(k) = getframe(h);
    k = k+1;
    clf()
    
end
close all

tStop = toc;

fprintf('time = %3.2f min\n',(tStop-tStart)/60)


writeVideo(writerObj, Frame);
close(writerObj);


%% Supplementary plot. Zone of inhibition for a given cell

% You can change this variable. It picks one cell in the lattice to analyze
% and show. The higher the number, the higher up lattice the cell is.
nCell = 25;

uI = cIJ(nCell,:);

uI = uI - min(uI);
uI = uI./max(uI);

figure(2)

axis([0 LX 0 LY])
axis equal


for nI = 1:nC
    
    patch(vX([c{nI}(1:end),c{nI}(1)]), vY([c{nI}(1:end),c{nI}(1)]),...
        [1 uI(nI) 0])
end

hold on
plot(latticeX(nCell),latticeY(nCell),'b*')

title('Interaction profile')
% The strength of inhibition from one cell is color coded
set(gcf, 'Units', 'Normalized', 'OuterPosition',...
    [0.1, 0.04, min(LX/LY,0.9), min(0.96,LY/LX)]);
axis([0 LX 0 LY])
axis equal


% 
% %% Supplementary video.  Shows the wave
% 
% 
% dt = 1;
% 
% for aL = 1 %[5, 7, 8, 10, 15, 20]
%     
%     k = 1;
%     
%     videoFps = 4;
%     %graymap = gray(256);
%     %colormap autumn
%     framePerImage = 1;
%     % movie for all images
%     writerObj = VideoWriter(...
%         sprintf('output/wave_%02d_fps_%02d_cell_T_%01.1f_a_%02.2f.mp4',...
%         videoFps,nN,aT,aL),'MPEG-4');
%     
%     writerObj.FrameRate = videoFps;
%     
%     open(writerObj);
%     
%     clear('Frame')
%     
%     for tI = 0:dt:2*T
%         
%         hF = figure(3);
%         set(gcf, 'Units', 'Normalized', 'OuterPosition',...
%             [0.1, 0.04, min(LX/LY,0.9), min(0.96,LY/LX)]);
%         
%         for nI = 1:nC
%             patch(vX([c{nI}(1:end),c{nI}(1)]), vY([c{nI}(1:end),c{nI}(1)]),...
%                 [1 activation_wave(tI,T,latticeY(nI),LY,aL,lD) 0])
%         end
%         drawnow;
%         %pause(0.3)
%         axis([0 LX 0 LY])
%         axis equal
%         
%         Frame(k) = getframe(hF);
%         k = k+1;
%         
%         fprintf('t = %02.2f\n',tI)
%         
%         clf();
%     end
%     
%     writeVideo(writerObj, Frame);
%     
%     close(writerObj);
%     
%     close(hF)
%     
% end

