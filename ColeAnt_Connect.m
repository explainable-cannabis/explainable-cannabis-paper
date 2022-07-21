%%script for computing average connectivity from a CIFTI file

%this script will
%1) read in the parcellation file from the Cole-Antecevic
%Cortical-Subcortical Networks
%2) read in the network assignment cifti for the Cole-Antecevic networks
%3) read in a cifti file for an individual
%4) average the timeseries within each of 718 parcels for that subject
%5) for each network, extract the nodal timeseries that are part of the
%network
%6) compute ROI-to-ROI connectivity using pearson correlations within each
%network
%7) z-transform and average the correlations within each network
%8) return 12 values for each subject corresponding to the average
%connectivity within each network
%
%Let's start!
%
%example use case:
%
%fmri_parc =('CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
%fmri_in =('rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii')
%net_assign =('cortex_subcortex_parcel_network_assignments.mat');
%conns =
%ELR_ColeAnt_Connect(path_to_functional,path_to_parcels,path_to_networks);

function Net_Corrs = ELR_ColeAnt_Connect(fmri_in, fmri_parc, net_assign)

%read in files, steps 1,2,3 above
fmri_in_parc=ft_read_cifti(fmri_parc);
fmri_func=ft_read_cifti(fmri_in);
netassignments=load(net_assign);

%find the weird nan-padding ft does
nan_pads = find(~isnan(fmri_in_parc.x1));

%remove the nan-padding
parcs = fmri_in_parc.x1(nan_pads);
fmri_timeseries = fmri_func.dtseries(nan_pads,:);

%average timeseries within each ROI, step 4
ROI_timeseries = [];
for i = unique(parcs)'
    temp_time = [];
    for j = 1:length(parcs)
        if parcs(j) == i
            temp_time = [temp_time; fmri_timeseries(j,:)];
        end
    end
    if length(find(parcs==i))>1
        ROI_timeseries = [ROI_timeseries; mean(temp_time)];
    else
        ROI_timeseries = [ROI_timeseries; temp_time];
    end   
end

%organize ROIs into 12 networks, step 5
%Compute Pearson correlations, step 6
Net_Corrs = [];
for i = 1:12
    temp_corr = [];
    for j = 1:718
        if netassignments(j)==i
            temp_corr = [temp_corr; ROI_timeseries(j,:)];
        end
    end
    [a,~]=size(temp_corr);
    corr_coeffs = [];
    for p = 1:a
        correl = corr(temp_corr(p,:)',temp_corr(p+1:end,:)');
        corr_coeffs = [corr_coeffs, correl];
    end
    %z-transform and average within network, step 7
    Net_Corrs = [Net_Corrs mean(atanh(corr_coeffs))];
end
%returns 12 correlation values!

end