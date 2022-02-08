clearvars; close all; 

%% process from raw data
strDir = 'D:\Raw data\USC System\20220203 Focus profile\'; 
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

keyboard; 


for nSet = 1 : 20
    % data directories     
%     strPrefix = sprintf('20220203_FocusDepth_%s_', nSet); 
    listFiles = dir(sprintf('%s20220203_FocusDepth_%d_*.dat', strDir, nSet)); 
    
    for nFile = 1 : 5
        strFile = sprintf('%s%s', strDir, listFiles(nFile).name); 



    end

    keyboard; 



end