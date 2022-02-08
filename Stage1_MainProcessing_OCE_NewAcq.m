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

for nSet = 22 % 21 : 50  %length(listRoot)
    
    
    % data directories
    
    strDir = ['C:\Users\NOIR_\Desktop\Jason\20220203 PDMS phantom test\Phantom_test_', num2str(nSet), '\']; 
    % strSaveMat = ['D:\Processed data\OCE and Scaffold\20220127 PDMS Phantom\', str];
    % strSaveFig = ['D:\Processed data\OCE and Scaffold\20220127 PDMS Phantom\', str, 'figures\'];
    % mkdir(strSaveMat); mkdir(strSaveFig); 
    
    listFolder = dir(sprintf('%s*.dat', strDir)); 
    
    
    
    % load calibration and dispersion
    strCalibration = 'D:\Calibration and Dispersion files\USC\20220127_USC_calibration.dat';
    [~, pdK, pnIndex] = readCalibration(strCalibration);
    pdK = pdK(1, :)';
    pnIndex = pnIndex(1, :)';
    
    strDispersion = 'D:\Calibration and Dispersion files\USC\20220127_USC_dispersion.dat'; 
    [pdDispReal, pdDispImag] = readDispersion(strDispersion); 
    
    
    % load header from the first data file
    strFile = fullfile(listFolder(1).folder, listFolder(1).name); 
    try
        cellArrays = readHeader(strFile); 
    catch
        cellArrays = readHeader2(strFile); 
    end
    
    nNumberLines = cellArrays{2,3};
    nLineLength = cellArrays{2,4};
    
    % processing data parameters
    dSurfThreshold = 80;    % surface calculation
    dK = 2 * pi / 560e-9; 
    dContrastMin = (1e6) * sqrt((4 / dK^2) * (1 - 10^(-0.5 / 10)));  % unit: um
    dContrastMax = (1e6) * sqrt((4 / dK^2) * (1 - 10^(-10 / 10)));
    
    
    % controling parameters
    bMATLABInterp = true; 
    bRef = 0;                % reference file: 1; no reference file: 0 
    pnNumbers = 51; %1 : 100; % numel(listData);
    
    %% calculate mask
    
    nLeft = 256;
    nRight = nLineLength - nLeft;
    nRound = 32;
    
    pdMask = calculateMask(nLineLength, nLeft, nRight, nRound);
    
    
    %% read reference spectrum / calculate reference spectrum
    if bRef         
        disp('reading reference...');
        listRef = dir(sprintf('%s*Ref*.dat', strDir));
        strReference = sprintf('%s%s', strDir, listRef(2).name);
        [pdIMAQ, ~] = readData(strReference, cellArrays); 
        pdRef = mean(pdIMAQ, 2);
    
%         % calculate noise level
%         pdRefSpectrum = pdIMAQ - repmat(pdRef, [1, size(pdIMAQ, 2)]); 
%         clear pdIMAQ; 
%         
%         pdRefCalibrated = applyCalibration(pdRefSpectrum, pdK, pnIndex, bMATLABInterp);
%         pcdRefCorrected = applyDispersion(pdRefCalibrated, pdDispReal, pdDispImag);
%         pcdRefDepthProfiles = getComplexDepthProfile(pcdRefCorrected, pdMask);
%         pcdRefDepthProfiles(round(nLineLength/2) + 1:end, : ) = []; 
%         clear pdRefCalibrated pcdRefCorrected; 
% 
%         pdNoise = mean(abs(pcdRefDepthProfiles) .^ 2, 2);
% 
%         clear pcdRefDepthProfiles;
    
    else 
        disp('calculating reference...');
        pdRef = -1;
        for nNumber = pnNumbers(1)+ 20 : -1 : pnNumbers(1)+2
            strFile = [listFolder(nNumber).folder, '\', listFolder(nNumber).name];
            disp(strFile); 
    
            [pdIMAQ, ~] = readData(strFile, cellArrays); 
            clear strFile; 
    
            if pdRef == -1
                pdRef = mean(pdIMAQ, 2);
            else
                pdRef = 0.9*pdRef + 0.1*mean(pdIMAQ, 2);
            end 
            clear pdIMAQ;
            
        end
        clear nNumber; 
    %     % noise level
    %     strFile = sprintf('%s%s', strDir, strReferenceFile);
    %     [pdIMAQ, ~] = readDataFile(strFile, cellArrays); 
    %     % calculate noise level
    %     pdSpectrum = pdIMAQ - repmat(pdRef, [1, size(pdIMAQ, 2)]); 
    %     pdSpectrum = correctSpectrum(pdSpectrum, cellArrays);
    %     clear pdIMAQ; 
    %     
    %     nMidLength = size(pdSpectrum, 1) / 2 + 1;
    %     pdFFT = fft(pdSpectrum);
    %     pdZPFFT = zeros([nZPPaddingFactor*size(pdSpectrum, 1), size(pdSpectrum, 2)]);
    %     pdFFT(nMidLength, :) = 0.5 * pdFFT(nMidLength, :);
    %     pdZPFFT(1:nMidLength, :) = pdFFT(1:nMidLength, :);
    %     pdZPFFT(end-nMidLength+2:end, :) = pdFFT(end-nMidLength+2:end, :);
    %     pdZPSpectrum = real(ifft(pdZPFFT)) * nZPPaddingFactor;
    %     clear nMidLength pdFFT pdZPFFT; 
    % 
    %     pnIndex = 2 * (0 : size(pdSpectrum, 1)-1);
    %     pdInterpolated = interp1(pdCalibrationZPPoint, pdZPSpectrum, pnIndex);
    %     pdCorrected = pdInterpolated .* repmat(pdDispersion', [1, size(pdInterpolated, 2)]);
    % 
    %     pdDepthProfiles = fft(pdCorrected);
    %     clear pdSpectrum pdZPSpectrum pnIndex pdInterpolated pdCorrected;
    %     pdDepthProfiles(size(pdDepthProfiles, 1)/2:end, :) = [];
    % 
    % %     pdRefDepthProfiles = pdDepthProfiles; 
    %     pdNoise = mean(abs(pdDepthProfiles), 2);
    end
    
    %% read actual data file
    disp('processing data...');
    
    for nNumber = pnNumbers
    
        % read IMAQ
        strFile = fullfile(listFolder(nNumber).folder, listFolder(nNumber).name); 
        fprintf('%d/%d: %s \n', nNumber, length(pnNumbers), strFile);
        [pdIMAQ, ~] = readData(strFile, cellArrays);
        
        %% processing intensity
        % subtract reference    
        pdIMAQ = pdIMAQ - repmat(pdRef,[1, nNumberLines]);
        
        % apply calibration   
        pdIMAQCalibrated = applyCalibration(pdIMAQ, pdK, pnIndex, bMATLABInterp);
        % clear pdIMAQ pdK pnIndex
        
        % apply dispersion
        pcdIMAQCorrected = applyDispersion(pdIMAQCalibrated, pdDispReal, pdDispImag);
        
        % get complex depth profile        
        pcdDepthProfiles = getComplexDepthProfile(pcdIMAQCorrected, pdMask);
        pcdDepthProfiles(round(nLineLength/2) + 1:end, : ) = []; 

        % intensity image        
        pdI = abs(pcdDepthProfiles) .^ 2; 
        pddB = 10* log10(pdI); 
    
        %% processing: surface
        pddBFilt = imfilter(pddB, ones([4, 16])/(4*16), 'replicate');
        pnSurface = zeros(1, size(pdI, 2)); 
    %     [pnX, pnY] = find(pdSmoothMu(101:end, :) > dSurfThreshold);
    %     pnSurfaceRaw(flipud(pnY)) = flipud(pnX)+100;
    %     pnSurface = round(smooth(pnSurfaceRaw, 20));    
    %     clear pnX pnY pnSurfaceRaw;
    
        [pnX, pnY] = find(pddBFilt(51:end, :) > dSurfThreshold);
        pnSurface(flipud(pnY)) = flipud(pnX) + 50;
    %     pnSurface = round(smooth(pnSurfaceRaw, 20)); 
        clear pnX pnY;
    
        %% processing: OCE fringe washout subtraction
        pddBOdd = imfilter(pddB(:, 1:2:end), ones([8, 8])/(8*8), 'replicate');
        pddBEven = imfilter(pddB(:, 2:2:end), ones([8, 8])/(8*8), 'replicate'); 
    
        pddBDiff = abs(pddBEven - pddBOdd);    
        pdDeltaZ = (1e6) * sqrt((4 / dK^2) * (1 - 10.^(-pddBDiff / 10)));

%         figure, [~, nROI] = imcrop(rescale(pdDeltaZ(1:1024, :))); 
%         pnX = round(nROI(1)) : round(nROI(1) + nROI(3)); 
%         pnY = round(nROI(2)) : round(nROI(2) + nROI(4)); 
% 
%         ddBDiff = mean(mean(pddBDiff(pnY, pnX))); 
%         dDeltaZ = mean(mean(pdDeltaZ(pnY, pnX))); 


        %% processing: attenuation
    
        %% show frame results and save data
        [~, strName, ~] = fileparts(listFolder(nNumber).name);
    
        figure(fA); 
        imagesc(pddB(1:1024, :), [55, 95]), colormap(1-gray), colorbar; 
%         hold on; plot(pnSurface, 'y', 'LineWidth', 2); hold off; 
        title(strName, 'interpreter', 'none'); 
%         saveas(gcf, sprintf('%sIntensity_%d_%s.png', strSaveFig, 100000 + nNumber, strName), 'png');

        figure(fD); 
        plot(mean(pddB(:, 501: 550), 2)); 
        xlim([1, 1024]); ylim([50, 100]);
    
        figure(fB); 
        imagesc(pddBDiff, [0, 10]); colormap(jet); colorbar; 
        hold on; plot(pnSurface(1:2:end), 'y', 'LineWidth', 2); hold off; 
%         hold on; rectangle('Position', nROI, 'EdgeColor', 'k', 'LineWidth', 2); hold off; 
%         title(sprintf('%s, dB diff = %0.2f', strName, ddBDiff), 'interpreter', 'none'); 
%         saveas(gcf, sprintf('%sOCE_%d_%s.png', strSaveFig, 100000 + nNumber, strName), 'png');
        
        figure(fC); 
        imagesc(pdDeltaZ, [dContrastMin, dContrastMax]); colormap(jet); colorbar; 
        hold on; plot(pnSurface(1:2:end), 'y', 'LineWidth', 2); hold off; 
%         hold on; rectangle('Position', nROI, 'EdgeColor', 'k', 'LineWidth', 2); hold off; 
%         title(sprintf('%s, del_z = %0.2f', strName, dDeltaZ), 'interpreter', 'none'); 
%         saveas(gcf, sprintf('%sDeltaZ_%d_%s.png', strSaveFig, 100000 + nNumber, strName), 'png');
    
%         save(sprintf('%s%d_%s.mat', strSaveMat, 100000 + nNumber, strName), 'pdDepthProfiles', 'pddB', 'pddBDiff', 'pnSurface', 'pdDeltaZ'); 
    
        
        drawnow;
    
    end

end
