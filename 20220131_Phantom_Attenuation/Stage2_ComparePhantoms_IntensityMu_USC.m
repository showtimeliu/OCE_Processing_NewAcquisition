clearvars; close all; 

% figure display parameters
scrsz = get(groot,'ScreenSize');
nY = 2;
nBottom = 50;
nTop = 90;
nX = 3;
nLeft = 10;
nRight = 10;
nHeight = scrsz(4)-nBottom;
nWidth = scrsz(3)-nLeft;
fA = figure('Position',[nLeft+0*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
fB = figure('Position',[nLeft+1*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
fC = figure('Position',[nLeft+2*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]); 
fD = figure('Position',[nLeft+0*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
fE = figure('Position',[nLeft+1*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
fF = figure('Position',[nLeft+2*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
clear scrsz nY nBottom nTop nX nLeft nRight nHeight nWidth;


strMatA = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_1\Phantom_test_1_100051.mat';
load(strMatA); 
pdIA = abs(pcdDepthProfiles .^ 2); 
pdMuA = pdMu; 
pdSmoothMuA = pdSmoothMu; 
pddBA = 10*log10(pdIA);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 

strMatB = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_5\Phantom_test_5_100051.mat';
load(strMatB); 
pdIB = abs(pcdDepthProfiles .^ 2); 
pdMuB = pdMu; 
pdSmoothMuB = pdSmoothMu; 
pddBB = 10*log10(pdIB);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 

strMatC = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_4\Phantom_test_4_100051.mat';
load(strMatC); 
pdIC = abs(pcdDepthProfiles .^ 2); 
pdMuC = pdMu; 
pdSmoothMuC = pdSmoothMu; 
pddBC = 10*log10(pdIC);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 

strMatD = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_6\Phantom_test_6_100051.mat';
load(strMatD); 
pdID = abs(pcdDepthProfiles .^ 2); 
pdMuD = pdMu; 
pdSmoothMuD = pdSmoothMu; 
pddBD = 10*log10(pdID);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 


strMatE = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_7\Phantom_test_7_100051.mat';
load(strMatE); 
pdIE = abs(pcdDepthProfiles .^ 2); 
pdMuE = pdMu; 
pdSmoothMuE = pdSmoothMu; 
pddBE = 10*log10(pdIE);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 


strMatF = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_8\Phantom_test_8_100051.mat';
load(strMatF); 
pdIF = abs(pcdDepthProfiles .^ 2); 
pdMuF = pdMu; 
pdSmoothMuF = pdSmoothMu; 
pddBF = 10*log10(pdIF);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 

strMatG = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Cartilage_test_1\Cartilage_test_1_100051.mat';
load(strMatG); 
pdIG = abs(pcdDepthProfiles .^ 2); 
pdMuG = pdMu; 
pdSmoothMuG = pdSmoothMu; 
pddBG = 10*log10(pdIG);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 

%% 
% A
nRangeA = [500, 600]; 
pddBProfileA = mean(pddBA(1:1024, nRangeA(1):nRangeA(2)), 2); 
pdXA = (((1:1024) - 255) * (2.0/1024))'; 
pdMuProfileA = mean(pdSmoothMuA(1:1024, nRangeA(1):nRangeA(2)), 2);

% B
nRangeB = [700, 800]; 
pddBProfileB = mean(pddBB(1:1024, nRangeB(1):nRangeB(2)), 2); 
pdXB = (((1:1024) - 198) * (2.0/1024))'; 
pdMuProfileB = mean(pdSmoothMuB(1:1024, nRangeB(1):nRangeB(2)), 2);

% C
nRangeC = [650, 750]; 
pddBProfileC = mean(pddBC(1:1024, nRangeC(1):nRangeC(2)), 2); 
pdXC = (((1:1024) - 181) * (2.0/1024))'; 
pdMuProfileC = mean(pdSmoothMuC(1:1024, nRangeC(1):nRangeC(2)), 2);

% D
nRangeD = [500, 600]; 
pddBProfileD = mean(pddBD(1:1024, nRangeD(1):nRangeD(2)), 2); 
pdXD = (((1:1024) - 269) * (2.0/1024))'; 
pdMuProfileD = mean(pdSmoothMuD(1:1024, nRangeD(1):nRangeD(2)), 2);

% E
nRangeE = [500, 600]; 
pddBProfileE = mean(pddBE(1:1024, nRangeE(1):nRangeE(2)), 2); 
pdXE = (((1:1024) - 207) * (2.0/1024))'; 
pdMuProfileE = mean(pdSmoothMuE(1:1024, nRangeE(1):nRangeE(2)), 2);

% F
nRangeF = [500, 600]; 
pddBProfileF = mean(pddBF(1:1024, nRangeF(1):nRangeF(2)), 2); 
pdXF = (((1:1024) - 207) * (2.0/1024))'; 
pdMuProfileF = mean(pdSmoothMuF(1:1024, nRangeF(1):nRangeF(2)), 2);

% G
nRangeG = [500, 600]; 
pddBProfileG = mean(pddBG(1:1024, nRangeG(1):nRangeG(2)), 2); 
pdXG = (((1:1024) - 496) * (2.0/1024))'; 
pdMuProfileG = mean(pdSmoothMuG(1:1024, nRangeG(1):nRangeG(2)), 2);


figure(fA); imagesc(pddBG, [50, 100]); colormap(1-gray); colorbar; 
ylim([1, 1024]);
title('Phantom F')

figure(fD); imagesc(pdMuG, [0, 5]); colormap(gray); colorbar; 
ylim([1, 1024]);
title('Phantom F')

figure(fE), plot(pddBProfileG);  

% % compare heating vs no heating curing 
% figure(fB); 
% pA = plot(pdXA, pddBProfileA, 'k', 'LineWidth', 2); 
% hold on; pB = plot(pdXB, pddBProfileB); 
% hold on; pC = plot(pdXC, pddBProfileC); 
% hold on; pD = plot(pdXD, pddBProfileD); 
% hold off; 
% xlim([-0.5, 1.5]); ylim([[50, 100]]); ylabel('intensity, dB'); xlabel('depth, mm');
% legend([pA pB pC pD], {'PDMS (1:10) no heat', 'PDMS (1:10) heat curing', 'PDMS (1:5) heat curing', 'PDMS (1:20) heat curing'}); 
% title('comparison of intensity profiles'); 
% 
% figure(fC); 
% pA = plot(pdXA, pdMuProfileA, 'k', 'LineWidth', 2);
% hold on; pB = plot(pdXB, pdMuProfileB);
% hold on; pC = plot(pdXC, pdMuProfileC);
% hold on; pD = plot(pdXD, pdMuProfileD);
% hold off; 
% xlim([-0.5, 1.5]); ylim([[0, 10]]); ylabel('attenuation'); xlabel('depth, mm');
% legend([pA pB pC pD], {'PDMS (1:10) no heat', 'PDMS (1:10) heat curing', 'PDMS (1:5) heat curing', 'PDMS (1:20) heat curing'}, ...
%     'Location', 'northwest'); 
% title('comparison of atteunation profiles'); 

% % compare PDMS vs silicone 
% figure(fB); 
% pA = plot(pdXB, pddBProfileB, 'b', 'LineWidth', 2); 
% hold on; pB = plot(pdXE, pddBProfileE); 
% hold on; pC = plot(pdXF, pddBProfileF); 
% hold off; 
% xlim([-0.5, 1.5]); ylim([[50, 100]]); ylabel('intensity, dB'); xlabel('depth, mm');
% legend([pA pB pC], {'PDMS (1:10) heat curing', 'silicone (hard)', 'silicone (soft)'}); 
% title('comparison of intensity profiles'); 
% 
% figure(fC); 
% pA = plot(pdXB, pdMuProfileB, 'b', 'LineWidth', 2);
% hold on; pB = plot(pdXE, pdMuProfileE);
% hold on; pC = plot(pdXF, pdMuProfileF);
% hold off; 
% xlim([-0.5, 1.5]); ylim([[0, 10]]); ylabel('attenuation'); xlabel('depth, mm');
% legend([pA pB pC], {'PDMS (1:10) heat curing', 'silicone (hard)', 'silicone (soft)'}, ...
%     'Location', 'northwest'); 
% title('comparison of atteunation profiles'); 


% compare PDMS vs silicone vs cartilage 
figure(fB); 
pA = plot(pdXB, pddBProfileB); 
hold on; pB = plot(pdXE, pddBProfileE); 
hold on; pC = plot(pdXG, pddBProfileG, 'LineWidth', 2); 
hold off; 
xlim([-0.5, 1.5]); ylim([[50, 100]]); ylabel('intensity, dB'); xlabel('depth, mm');
legend([pA pB pC], {'PDMS (1:10) heat curing', 'silicone (hard)', 'bovine cartilage (control)'}); 
title('comparison of intensity profiles'); 

figure(fC); 
pA = plot(pdXB, pdMuProfileB);
hold on; pB = plot(pdXE, pdMuProfileE);
hold on; pC = plot(pdXG, pdMuProfileG, 'LineWidth', 2);
hold off; 
xlim([-0.5, 1.5]); ylim([[0, 10]]); ylabel('attenuation'); xlabel('depth, mm');
legend([pA pB pC], {'PDMS (1:10) heat curing', 'silicone (hard)', 'bovine cartilage (control)'}, ...
    'Location', 'northwest'); 
title('comparison of atteunation profiles'); 

