clear;

d0=pwd;

subjects = {'Br2934', 'Br5347', 'Br5394', 'Br5766', 'Br5785','Br9413'};
%subjects = {'Br2934'};

plotDir = fullfile(d0,'../../plots/08_validation_oligos/');
dataDir = fullfile(d0,'../../processed_data/08_validation_oligos/');

layType = 'equivol';

% set up for visualisations

% % deinfe corresponding colomaps (like FSL)
% cmapR = [linspace(101/256,255,256)',zeros(256,2)];
% cmapBLB = [zeros(256,1),linspace(0,1,256)',ones(256,1)];
% cmapRY = [ones(256,1),linspace(0,1,256)',zeros(256,1)];
% cmapG = [zeros(256,1),linspace(100/256,1,256)',zeros(256,1)];
% cmaps = {cmapR,cmapB,cmapY,cmapG,cmapG,cmapG,gray};

% define colourblind friendly colomaps
cmapO = [linspace(0,238/255,256)',linspace(0,119/255,256)',linspace(0,51/255,256)'];
cmapLGrey = repmat(linspace(0,0.5,256)',1,3);
cmapB = [linspace(0,68/255,256)',linspace(0,119/255,256)',linspace(0,170/255,256)'];
cmapY = [linspace(0,238/255,256)',linspace(0,221/255,256)',linspace(0,136/255,256)'];
cmapG = [zeros(256,1),linspace(100/256,1,256)',zeros(256,1)];
cmapLay = [238 102 119; 204 187 68; 34 136 51; 102 204 238]./255;
cmaps = {cmapO,cmapB,cmapY,cmapG,cmapG,cmapLay,gray};

% define ranges
clims = {[1,5],[1,5],[1,5],[1,2],[1,4],[1,4],[1,2]};

% % define groups for each plot (7 is background)
% plotGroups = {[7 1 2 3],...
%     [7 1 2 3 4],...
%     [7 1 2 3 5],...
%     [7 1 2 3 6]};

% maybe just show original and final
plotGroups = {[7 1 2 3],[7 1 2 3 6]};

f1 = figure('Color','k','Position',get(0,'ScreenSize'));
t = tiledlayout(length(plotGroups),length(subjects),...
        'Padding','tight','TileSpacing','tight');

for S = 1:length(subjects)

    % First create a masking visualisation
    % Raw data (DAPI DAPI positive; MBP Oligo; HPCAL1_L2_neuron)
    % Boundaries
    % Layers
    % Layers masked

    thisSub = subjects{S};

    disp(['+++++ ' thisSub ' +++++']);

    subDir = fullfile(plotDir, thisSub);

    % load in images
    visNiiNames ={...
        [thisSub '_array_DAPI_DAPI_Positive_128x128binned_for_masking.nii'],...
        [thisSub '_array_MBP_Oligo_128x128binned_for_masking.nii'],...
        [thisSub '_array_HPCAL1_L2_neuron_128x128binned_for_masking.nii'],...
        [thisSub '_128x128binned_manual_boundaries.nii'],...
        [thisSub '_128x128binned_laynii_layers_' layType '.nii'],...
        [thisSub '_128x128binned_laynii_layers_' layType '_masked.nii']
        };
    for ii = 1:length(visNiiNames)
        visNii = nifti(fullfile(subDir, visNiiNames{ii}));
        visMtxs{ii} = visNii.dat(:,:,:);
    end
    % define a background
    visMtxs{7} = ones(size(visMtxs{1}));
    clear visNii

    for p = 1:length(plotGroups)

        axBase = nexttile(t,S + length(subjects)*(p-1));
        axis off;

        hold on;
        
        if p == 1
            title(thisSub,'Color','w','FontSize',14,'HorizontalAlignment','center')
        end

        pgroup = plotGroups{p};

        for ii = 1:length(pgroup)

            idx = pgroup(ii);

            ax = axes('Position',axBase.Position);
            imagesc(visMtxs{idx},'Parent',ax,clims{idx});
            colormap(ax,cmaps{idx});
            
            if ismember(ii,[5,6])
                alphaDat = double(visMtxs{idx} ~= 0)*0.8;
            else 
                alphaDat = double(visMtxs{idx} ~= 0);
            end

            set(get(ax,'Children'),'AlphaData',alphaDat);

            axis off;
            set(ax,'Color','none')
            axis image

        end

    end

    % Analyses
    % Oligo count in different layers (normalised by size of layer)
    % Oligo cytoplasm intensity (normalised by overall intensity within gm mask)

    % load in matrices to be analysed
    dapiCount = rgb2gray(imread(...
        fullfile(subDir,[thisSub '_array_1x1binned_DAPI (DAPI) Positive.png'])));
    oligoCount = rgb2gray(imread(...
        fullfile(subDir,[thisSub '_array_1x1binned_consensus_MBP.png'])));
    %intensity calculations will need to be done using a different method
    % mbpCyto = rgb2gray(imread(...
    %     fullfile(subDir,[thisSub '_array_1x1binned_MBP (Opal 520) Cytoplasm Intensity.png'])));
    % mbpNucl = rgb2gray(imread(...
    %     fullfile(subDir,[thisSub '_array_1x1binned_MBP (Opal 520) Nucleus Intensity.png'])));

    % upscale lay128bin
    lay128bin = uint8(visMtxs{6});
    lay1bin = imresize(lay128bin,[size(oligoCount,1),size(oligoCount,2)],'nearest');
    % save a .txt file for the 1x1 binned SpD images 
    writematrix(lay1bin, fullfile(subDir, [thisSub '_1x1binned_laynii_layers_' layType '_masked.txt']));

    % load in and upscale white matter
    wm128bin = nifti(fullfile(subDir,[thisSub '_128x128binned_manual_white_matter_masked.nii']));
    wm128bin = uint8(wm128bin.dat(:,:,:));
    wm1bin = imresize(wm128bin,[size(oligoCount,1),size(oligoCount,2)],'nearest');

    % example oligo counts
    for SpD = 1:4
        % raw oligo count 
        SpDOligCountRaw(SpD,S) =...
            bwconncomp(oligoCount & lay1bin == SpD, 4).NumObjects;
        % proportion of oligos
        SpDOligCountNorm(SpD,S) =...
            100*SpDOligCountRaw(SpD,S)/...
            (bwconncomp(dapiCount & lay1bin == SpD, 4).NumObjects);
        % % MBP intensity in cyto
        % SpDMBPCytoRaw(SpD,S) =...
        %     median(mbpCyto(mbpCyto(:) ~= 0 & lay1bin(:) == SpD));
        % SpDMBPCytoNorm(SpD,S) =...
        %     double(SpDMBPCytoRaw(SpD,S))/...
        %     double(median(mbpCyto(mbpCyto(:) ~= 0 & wm1bin(:) ~= 0)));
        % % MBP intensity in nucl
        % SpDMBPNuclRaw(SpD,S) =...
        %     median(mbpNucl(mbpNucl(:) ~= 0 & lay1bin(:) == SpD));
        % SpDMBPNuclNorm(SpD,S) =...
        %     double(SpDMBPNuclRaw(SpD,S))/...
        %     double(median(mbpNucl(mbpNucl(:) ~= 0 & wm1bin(:) ~= 0)));
        fprintf('> SpD%d\n', SpD);
        fprintf('oligo count: %d\n', SpDOligCountRaw(SpD,S));
        fprintf('oligo proportion: %.4g\n', SpDOligCountNorm(SpD,S));
        % fprintf('MBP cyto raw intensity: %.4g\n', SpDMBPCytoRaw(SpD,S));
        % fprintf('MBP cyto norm intensity: %.4g\n', SpDMBPCytoNorm(SpD,S));
        % fprintf('MBP nucl raw intensity: %.4g\n', SpDMBPNuclRaw(SpD,S));
        % fprintf('MBP nucl norm intensity: %.4g\n', SpDMBPNuclNorm(SpD,S));

    end

end

set(f1, 'InvertHardCopy', 'off');
set (f1,'Color','k');
saveas(f1,fullfile(plotDir,['SpD_segmentations_' layType '.png']));

f2 = figure();
subplot(1,2,1)
plot(SpDOligCountRaw','x', 'MarkerSize',10)
xlim([0,7])
xticks(1:6)
xticklabels(subjects)
title('Raw olig. count')
subplot(1,2,2)
plot(SpDOligCountNorm','x', 'MarkerSize',10)
xlim([0,7])
xticks(1:6)
xticklabels(subjects)
legend({'SpD1','SpD2','SpD3','SpD4'},'Location','best')
title('Olig. proportion (%)')
saveas(f2,fullfile(plotDir,['SpD_oligos_by_subject_' layType '.png']));

% f3 = figure();
% subplot(1,2,1)
% plot(SpDMBPCytoRaw','x','MarkerSize',10)
% xlim([0,7])
% xticks(1:6)
% xticklabels(subjects)
% title('Raw MBP cyto int.')
% subplot(1,2,2)
% plot(SpDMBPCytoNorm','x','MarkerSize',10)
% xlim([0,7])
% xticks(1:6)
% xticklabels(subjects)
% title('Norm. MBP cyto int.')
% legend({'SpD1','SpD2','SpD3','SpD4'},'Location','best')
% saveas(f3,fullfile(plotDir,['SpD_MBP_cyto_by_subject_' layType '.png']));
% 
% f4 = figure();
% subplot(1,2,1)
% plot(SpDMBPNuclRaw','x','MarkerSize',10)
% xlim([0,7])
% xticks(1:6)
% xticklabels(subjects)
% title('Raw MBP nucl int.')
% subplot(1,2,2)
% plot(SpDMBPNuclNorm','x','MarkerSize',10)
% xlim([0,7])
% xticks(1:6)
% xticklabels(subjects)
% title('Norm. MBP nucl int.')
% legend({'SpD1','SpD2','SpD3','SpD4'},'Location','best')
% saveas(f4,fullfile(plotDir,['SpD_MBP_nucl_by_subject_' layType '.png']));

SpDRawOligoData = [table(subjects','VariableNames',{'Subject'}),...
    array2table(SpDOligCountRaw','VariableNames',...
    {'OligCountRaw_SpD1',...
    'OligCountRaw_SpD2',...
    'OligCountRaw_SpD3',...
    'OligCountRaw_SpD4'})];
writetable(SpDRawOligoData,...
    fullfile(dataDir,['SpD_oligo_count_raw_' layType '.csv']));

SpDNormOligoData = [table(subjects','VariableNames',{'Subject'}),...
    array2table(SpDOligCountNorm','VariableNames',...
    {'OligCountNorm_SpD1',...
    'OligCountNorm_SpD2',...
    'OligCountNorm_SpD3',...
    'OligCountNorm_SpD4'})];
writetable(SpDNormOligoData,...
    fullfile(dataDir,['SpD_oligo_count_norm_' layType '.csv']));
