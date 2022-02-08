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

for nSet = 7 % 21 : 50  %length(listRoot)
    
    
    % data directories
    str = ['20220203 PDMS phantom attenuation\Phantom_test_', num2str(nSet), '\']; 
    strDir = ['C:\Users\NOIR_\Desktop\Jason\', str]; 
    strSaveMat = ['D:\Processed data\OCE and Scaffold\', str];
    strSaveFig = ['D:\Processed data\OCE and Scaffold\', str, 'figures\'];
%     mkdir(strSaveMat); mkdir(strSaveFig); 
    
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
    dPixelDepth = 1; % 4.0 / (nLineLength/2); 
    
    
    % controling parameters
    bMATLABInterp = true; 
    bRef = 1;                % reference file: 1; no reference file: 0 
    pnNumbers = 50; %1 : 100; % numel(listData);
    
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
    
        % calculate noise level
        pdRefSpectrum = pdIMAQ - repmat(pdRef, [1, size(pdIMAQ, 2)]); 
        clear pdIMAQ; 
        
        pdRefCalibrated = applyCalibration(pdRefSpectrum, pdK, pnIndex, bMATLABInterp);
        pcdRefCorrected = applyDispersion(pdRefCalibrated, pdDispReal, pdDispImag);
        pcdRefDepthProfiles = getComplexDepthProfile(pcdRefCorrected, pdMask);
        pcdRefDepthProfiles(round(nLineLength/2) + 1:end, : ) = []; 
        clear pdRefCalibrated pcdRefCorrected; 

        pdNoise = mean(abs(pcdRefDepthProfiles) .^ 2, 2);

        clear pcdRefDepthProfiles;
    
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
        
        pdDepth = mean(pddB(:, 501: 550), 2); 
        pdX = (((1:2048) - 199) * (4.0/2048))'; 

        pnInd = find(pdX > 0.15 & pdX < 0.8); 
        p = polyfit(pdX(pnInd), pdDepth(pnInd), 1); 
        pdFit = polyval(p, pdX(pnInd)); 

    
        %% processing: OCE fringe washout subtraction
        pddBOdd = imfilter(pddB(:, 1:2:end), ones([8, 8])/(8*8), 'replicate');
        pddBEven = imfilter(pddB(:, 2:2:end), ones([8, 8])/(8*8), 'replicate'); 
    
        pddBDiff = abs(pddBEven - pddBOdd);    
        pdDeltaZ = (1e6) * sqrt((4 / dK^2) * (1 - 10.^(-pddBDiff / 10)));


        %% processing: attenuation
        pdIDenoised = pdI - repmat(pdNoise, [1, size(pdI, 2)]); 
        pdIDenoised(pdIDenoised <= 0) = 10e-7; 

        % initialize mu and sum
        pdMu = 0 * pdI; 
        pdSum = zeros(1, size(pdI, 2)); 
        % calculate mu; 
        for nDepth = (size(pdI, 1)-1) : -1 : 1
            pdSum = pdSum + pdIDenoised(nDepth+1, :); 
            pdMu(nDepth, :) = (1 / (2 * dPixelDepth)) * log(1 + pdIDenoised(nDepth, :) ./ pdSum);         
        end
        clear pdIDenoised pdSum nDepth; 

%         pdMu2 = imfilter(pdMu, ones(16, 16)/(16*16), 'replicate');

        % calculate smooth mu
        % subtract noise
        pdIDenoised = medfilt2(pdI, [16, 32]) - repmat(pdNoise, [1, size(pdI, 2)]); 
        pdIDenoised(pdIDenoised <= 0) = 10e-7; 
        % initialize mu and sum
        pdSmoothMu = 0 * pdI; 
        pdSum = zeros(1, size(pdI, 2)); 
        % calculate mu; 
        for nDepth = (size(pdI, 1)-1) : -1 : 1
            pdSum = pdSum + pdIDenoised(nDepth+1, :); 
            pdSmoothMu(nDepth, :) = (1 / (2 * dPixelDepth)) * log(1 + pdIDenoised(nDepth, :) ./ pdSum);         
        end
        clear pdIDenoised pdSum nDepth;



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




    
        %% show frame results and save data
        [~, strName, ~] = fileparts(listFolder(nNumber).name);
    
        figure(fA); 
        imagesc(pddB, [55, 95]), colormap(1-gray), colorbar; 
%         hold on; plot(pnSurface, 'y', 'LineWidth', 2); hold off; 
        title(strName, 'interpreter', 'none'); 
%         saveas(gcf, sprintf('%sIntensity_%d_%s.png', strSaveFig, 100000 + nNumber, strName), 'png');

        figure(fD); 
        plot(pdX, pdDepth); 
        hold on; plot(pdX(pnInd), pdFit, 'b', 'LineWidth', 2); hold off; 
        xlim([-0.5, 1.5]); ylim([50, 100]);
        title(sprintf('slope: %0.4f', p(1)));

        figure(fE); 
        plot(10*log10(mean(pdMu(:, 501:550), 2))); 

        figure(fF); 
        plot(10*log10(mean(pdSmoothMu(:, 501:550), 2))); 

    
%         figure(fB); 
%         imagesc(pddBDiff(1:1024, :), [0, 10]); colormap(jet); colorbar; 
%         hold on; plot(pnSurface(1:2:end), 'y', 'LineWidth', 2); hold off; 
%         title(strName, 'interpreter', 'none'); 
% %         saveas(gcf, sprintf('%sOCE_%d_%s.png', strSaveFig, 100000 + nNumber, strName), 'png');
%         
%         figure(fC); 
%         imagesc(pdDeltaZ(1:1024, :), [dContrastMin, dContrastMax]); colormap(jet); colorbar; 
%         hold on; plot(pnSurface(1:2:end), 'y', 'LineWidth', 2); hold off; 
%         title(strName, 'interpreter', 'none'); 
%         saveas(gcf, sprintf('%sDeltaZ_%d_%s.png', strSaveFig, 100000 + nNumber, strName), 'png');

        figure(fB); 
        imagesc(10*log10(pdMu), [-40, -20]); colormap(gray); colorbar;
        title(sprintf('%s, mu', strName), 'interpreter', 'none'); 
    
        figure(fC); 
        imagesc(10*log10(pdSmoothMu), [-40, -20]); colormap(gray); colorbar; 
        title(sprintf('%s, smooth mu', strName), 'interpreter', 'none'); 
    
%         save(sprintf('%s%s.mat', strSaveMat, strName), 'pcdDepthProfiles', 'pddBDiff', 'pnSurface', 'pdMu', 'pdSmoothMu', 'pdNoise'); 
%         save(sprintf('%s%d_%s.mat', strSaveMat, 100000 + nNumber, strName), 'pdDepthProfiles', 'pddB', 'pddBDiff', 'pnSurface', 'pdDeltaZ'); 
    
        
        drawnow;
    
    end

end
