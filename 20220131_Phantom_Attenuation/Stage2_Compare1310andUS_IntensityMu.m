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

strMat1310 = 'D:\Processed data\OCE and Scaffold\20220125 PDMS phantom under 1310 system\PDMS_phantom_1_1310_water\20220125_PDMS_phantom_1_1310_water_100275.mat';
load(strMat1310); 

pdI_1310 = circshift(pdIntensity, [0, -59]); 
pdMu_1310 = circshift(pdMu, [0, -59]);
pdSmoothMu_1310 = circshift(pdSmoothMu, [0, -59]);
pddB_1310 = 10*log10(pdI_1310);
clear pdIntensity pcdParallelOdd pcdParallelEven pcdPerpendicularOdd pcdPerpendicularEven; 
clear pdMu pdSmoothMu; 

strMat560 = 'D:\Processed data\OCE and Scaffold\20220128 PDMS phantom test\Phantom_test_1\Phantom_test_1_100051.mat';
load(strMat560); 
pdI_560 = abs(pcdDepthProfiles .^ 2); 
pdMu_560 = pdMu; 
pdSmoothMu_560 = pdSmoothMu; 
pddB_560 = 10*log10(pdI_560);
clear pcdDepthProfiles pddBDiff pdMu pdNoise pdSmoothMu pnSurface; 


% depth profiles 
nRange1310 = [150, 200]; 
pddBProfile_1310 = mean(pddB_1310(:, nRange1310(1):nRange1310(2)), 2); 
pdX1 = (((1:512) - 124) * (2.0/512))'; 
pnInd1 = find(pdX1 > 0.25 & pdX1 < 1.0); 
p1 = polyfit(pdX1(pnInd1), pddBProfile_1310(pnInd1), 1); 
pdFit1 = polyval(p1, pdX1(pnInd1)); 

nRange560 = [500, 600]; 
pddBProfile_560 = mean(pddB_560(1:1024, nRange560(1):nRange560(2)), 2); 
pdX2 = (((1:1024) - 255) * (2.0/1024))'; 
pnInd2 = find(pdX2 > 0.25 & pdX2 < 1.0); 
p2 = polyfit(pdX2(pnInd2), pddBProfile_560(pnInd2), 1); 
pdFit2 = polyval(p2, pdX2(pnInd2)); 


% mu depth profiles 
pdMuProfile_1310 = mean(pdMu_1310(:, nRange1310(1):nRange1310(2)), 2);
pdMuProfile_560 = mean(pdMu_560(1:1024, nRange560(1):nRange560(2)), 2);





% image display
figure(fA); 
imagesc(pddB_1310, [50, 100]); colormap(1-gray); colorbar; 
hold on; plot([nRange1310(1), nRange1310(1)], [1, 1024], 'r--'); 
hold on; plot([nRange1310(2), nRange1310(2)], [1, 1024], 'r--');
hold off; 
title('1310 nm, phantom 1, intensity'); 

figure(fB); 
imagesc(pddB_560, [50, 100]); colormap(1-gray); colorbar; 
hold on; plot([nRange560(1), nRange560(1)], [1, 1024], 'r--'); 
hold on; plot([nRange560(2), nRange560(2)], [1, 1024], 'r--');
ylim([1, 1024]); 
title('560 nm, phantom 1, intensity'); 

figure(fC); 
pl1 = plot(pdX1, pddBProfile_1310); 
hold on; pl2 = plot(pdX2, pddBProfile_560);
hold on; plot(pdX1(pnInd1), pdFit1, 'b', 'LineWidth', 2); 
hold on; plot(pdX2(pnInd2), pdFit2, 'r', 'LineWidth', 2);
hold off; 
xlim([-0.5, 1.5]); ylim([[50, 100]]); ylabel('intensity, dB'); xlabel('depth, mm');
legend([pl1 pl2], {'1310 nm system', '560 nm system'}); 
title('comparison of intensity profiles'); 

figure(fD);
imagesc(pdMu_1310, [0, 5]); colormap(gray); colorbar; 
hold on; plot([nRange1310(1), nRange1310(1)], [1, 1024], 'r--'); 
hold on; plot([nRange1310(2), nRange1310(2)], [1, 1024], 'r--');
hold off; 
title('1310 nm, phantom 1, attenuation'); 

figure(fE);
imagesc(pdMu_560, [0, 5]); colormap(gray); colorbar;
hold on; plot([nRange560(1), nRange560(1)], [1, 1024], 'r--'); 
hold on; plot([nRange560(2), nRange560(2)], [1, 1024], 'r--');
ylim([1, 1024]); 
title('560 nm, phantom 1, attenuation'); 

figure(fF);
pl1 = plot(pdX1, pdMuProfile_1310); 
hold on; pl2 = plot(pdX2, pdMuProfile_560);
xlim([-0.5, 1.5]); ylim([[0, 10]]); ylabel('attenuation'); xlabel('depth, mm');
legend([pl1 pl2], {'1310 nm system', '560 nm system'}); 
title('comparison of atteunation profiles'); 

