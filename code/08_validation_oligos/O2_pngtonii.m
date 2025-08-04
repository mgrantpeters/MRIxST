addpath('../../../PM_validation/nii_tool/')

%subjects = {'Br2934', 'Br5347', 'Br5394', 'Br5766', 'Br5785', 'Br9413'};
subjects = {'Br9413'};

plotDir = fullfile('../../plots/08_validation_oligos/');

nii_tool('RGBstyle','fsl');

% load in png files masking
for ii = 1:length(subjects)
    
    thisSub = subjects{ii};
    disp(['>>> ' thisSub]);
        
    subPlotDir = fullfile(plotDir, thisSub);

    imInfo = dir(fullfile(subPlotDir,'*128x128binned*.png'));

    for jj = 1:length(imInfo)
        
        img = imread(fullfile(subPlotDir, imInfo(jj).name));
        disp(imInfo(jj).name)
        % convert to nii
        nii = nii_tool('init',img);
        nii.img = rgb2gray(nii.img);

        disp(['size: ' num2str(size(nii.img))]);
        
        newName = strrep(strrep(strrep(strrep(imInfo(jj).name, ' (', '_'),') ', '_'),' ', '_'), '.png','.nii');

        nii_tool('save',nii, fullfile(subPlotDir,newName));

    end
end