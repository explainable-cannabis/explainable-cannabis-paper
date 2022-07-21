%code to parcellate an fMRI file according to an atlas in nifti format
%In current form, both files must be in the same resolution

%argument 1: input niftii file
%argument 2: input atlas file
%argument 3: output filename
%argument 4: do PCA? if 1, then do PCA (otherwise just average)
%argument 5: if PCA, how much variance to keep? then []

function [ROIfunc,PercZeros] = ELR_ParcellateNifti(infile,atlas,outfile,doPCA,PCAvariance)

NiifMRI=niftiread(infile);
ROIs=niftiread(atlas);
if size(squeeze(NiifMRI(:,:,:,1))) ~= size(ROIs)
    error('fMRI and atlas are not the same resolution! Aborting.');
end
%loop through ROIs and assign voxels for averaging
ROI_list = unique(ROIs);
ROI_list = ROI_list(2:end);
NumROIs = numel(ROI_list);
[~,~,~,time]=size(NiifMRI);
ROIfunc = zeros(NumROIs,time);
PercZeros=zeros(NumROIs,1);
for k = 1:NumROIs
    %find linear indices corresponding to voxels in the ROI, then
    %convert to matrix indices
    [ROIindX,ROIindY,ROIindZ] = ind2sub(size(NiifMRI),find(ROIs==ROI_list(k)));
    %populate a temporary file with the timeseries from those
    %voxels
    ROItemp = zeros(numel(ROIindX),time);
    for j = 1:length(ROIindX)
        ROItemp(j,:) = NiifMRI(ROIindX(j),ROIindY(j),ROIindZ(j),:);
    end
    %find number of zeros (voxels in ROI but not in fMRI)
    PercZeros(k) = numel(find(mean(ROItemp,2)==0))/numel(ROIindX);
    %if a ROI is not in this subject skip it...
    if PercZeros(k) == 1
        ROIfunc(k,:) = zeros(1,time);
        break
    end
    ROItemp(mean(ROItemp,2)==0,:)=[];
    %average over the voxels OR do PCA and keep the first n components
    %that explain v% variance 
    if doPCA == 1
        ROItemp=ROItemp';
        warning off
        [~,ROIPCA,~,~,vars] = pca(ROItemp);
        warning on
        %find the necessary number of components for specified amount of
        %variance
        varsum=0;
        for h = 1:numel(vars)
            varsum = varsum+vars(h);
            if varsum >= PCAvariance
                varneeded = h;
                break
            end
        end
        ROIPCA = ROIPCA.*vars';
        ROIfunc(k,:) = mean(ROIPCA(:,1:varneeded),2)';
    else
        ROIfunc(k,:) = mean(ROItemp);
    end
end
%save this as a .mat output
save(outfile,'ROIfunc');
save([outfile '_NumZeros'],'PercZeros');
end
