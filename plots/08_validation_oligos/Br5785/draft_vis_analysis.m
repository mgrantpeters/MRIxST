
% First create a masking visualisation
% Raw data (DAPI DAPI positive; MBP Oligo; HPCAL1_L2_neuron)
% Boundaries
% Layers
% Layers masked

% load in images
visNiiNames = {'Br5785_array_DAPI_DAPI_Positive_128x128binned_for_masking.nii',...
    'Br5785_array_MBP_Oligo_128x128binned_for_masking.nii',...
    'Br5785_array_HPCAL1_L2_neuron_128x128binned_for_masking.nii',...
    'Br5785_128x128binned_manual_boundaries.nii',...
    'Br5785_128x128binned_laynii_layers_equivol.nii',...
    'Br5785_128x128binned_laynii_layers_equivol_masked.nii'
    };
for ii = 1:length(visNiiNames)
    visNii = nifti(visNiiNames{ii});
    visMtxs{ii} = visNii.dat(:,:,:);
end
% define a background
visMtxs{7} = ones(size(visMtxs{1}));
clear visNii

% deinfe corresponding colomaps (like FSL)
cmapR = [linspace(100/256,1,256)',zeros(256,2)];
cmapBLB = [zeros(256,1),linspace(0,1,256)',ones(256,1)];
cmapRY = [ones(256,1),linspace(0,1,256)',zeros(256,1)];
cmapG = [zeros(256,1),linspace(100/256,1,256)',zeros(256,1)];
cmaps = {cmapR,cmapBLB,cmapRY,cmapG,cmapG,cmapG,gray};

% define ranges
clims = {[1,5],[1,5],[1,5],[1,2],[1,4],[1,4],[1,2]};

% % define groups for each plot (7 is background)
% plotGroups = {[7 1 2 3],...
%     [7 1 2 3 4],...
%     [7 1 2 3 5],...
%     [7 1 2 3 6]};

% maybe just show original and final
plotGroups = {[7 1 2 3],[7 1 2 3 6]};

f1 = figure('Color','k');
t = tiledlayout(length(plotGroups),1,'Padding','tight','TileSpacing','tight');

for p = 1:length(plotGroups)
    
    axBase = nexttile(t);
    axis off;

    hold on;

    pgroup = plotGroups{p};
    
    for ii = 1:length(pgroup)
        
        idx = pgroup(ii);

        ax = axes('Position',axBase.Position);
        imagesc(visMtxs{idx},'Parent',ax,clims{idx});
        colormap(ax,cmaps{idx});

        alphaDat = double(visMtxs{idx} ~= 0);
        set(get(ax,'Children'),'AlphaData',alphaDat);

        axis off;
        set(ax,'Color','none')
        axis image
    
    end

end

% Analyses
% Oligo count in different layers (normalised by size of layer)
% Oligo cytoplasm intensity (normalised by overall intensity within gm mask)

% upscale lay128bin
lay128bin = uint8(visMtxs{6});
lay1bin = repelem(lay128bin,128,128);

% load in matrices to be analysed
oligoCount = rgb2gray(imread('Br5785_array_1x1binned_MBP_Oligo.png'));
mbpCyto = rgb2gray(imread('Br5785_array_1x1binned_MBP (Opal 520) Cytoplasm Intensity.png'));

% example oligo counts
for SpD = 1:4 
    SpDOligCountRaw(SpD) = sum(logical((oligoCount(lay1bin == SpD))));
    SpDOligCountPerPix(SpD) = SpDOligCountRaw(SpD)/sum(lay1bin(:) == SpD);
    SpDMBPCytoRaw(SpD) = mean(mbpCyto(lay1bin(:) == SpD));
    SpDMBPCytoNorm(SpD) = SpDMBPCytoRaw(SpD)/mean(mbpCyto(lay1bin(:) ~= 0));
    fprintf('> SpD%d\n', SpD);
    fprintf('oligo count: %d\n', SpDOligCountRaw(SpD));
    fprintf('oligos per pixel: %.4g\n', SpDOligCountPerPix(SpD));
    fprintf('MBP cyto raw intensity: %.4g\n', SpDMBPCytoRaw(SpD));
    fprintf('MBP cyto norm intensity: %.4g\n', SpDMBPCytoNorm(SpD));
end

